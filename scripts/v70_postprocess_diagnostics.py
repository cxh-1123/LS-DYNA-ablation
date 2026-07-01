"""
v70_postprocess_diagnostics.py
==============================

V7 diagnostics for existing LAMMPS TTM-MD Si ablation trajectories.

This script does not modify or rerun the simulation.  It reads the existing
V61 ttm/mod trajectory, builds physically motivated diagnostics, and produces
continuous x-z field plots suitable for deciding what the next simulation must
change.

Example:
    python scripts\\v70_postprocess_diagnostics.py --project-root D:\\cxh-daima\\LS-DYNA-ablation --case fluence_0p5_100ps
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import deque
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


KB_EV_K = 8.617333262145e-5
SI_MELT_K = 1687.0
SI_VAP_K = 3538.0
ANG_TO_NM = 0.1
BAR_TO_GPA = 1.0e-4


def iter_lammpstrj(path: Path):
    with path.open("r", encoding="utf-8", errors="replace") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(f.readline().strip())
            assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
            natoms = int(f.readline().strip())
            bounds_header = f.readline().strip()
            if not bounds_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Unexpected bounds header: {bounds_header}")
            bounds = []
            for _ in range(3):
                lo, hi, *_ = f.readline().split()
                bounds.append((float(lo), float(hi)))
            atom_header = f.readline().strip()
            if not atom_header.startswith("ITEM: ATOMS"):
                raise ValueError(f"Unexpected atom header: {atom_header}")
            cols = atom_header.split()[2:]
            data = np.loadtxt([f.readline() for _ in range(natoms)])
            if data.ndim == 1:
                data = data.reshape(1, -1)
            yield step, bounds, cols, data


def parse_thermo(log_path: Path) -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    if not log_path.exists():
        return rows
    header: list[str] | None = None
    with log_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "Step" and parts[1] == "Time":
                header = parts
                continue
            if header and len(parts) >= len(header):
                try:
                    vals = [float(x) for x in parts[: len(header)]]
                except ValueError:
                    header = None
                    continue
                row = dict(zip(header, vals))
                rows[int(row["Step"])] = row
    return rows


def read_te_file(path: Path) -> tuple[float, float, float, float]:
    vals = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                vals.append(float(parts[3]))
    arr = np.asarray(vals, dtype=float)
    if arr.size == 0:
        return math.nan, math.nan, math.nan, math.nan
    return (
        float(np.nanmin(arr)),
        float(np.nanmean(arr)),
        float(np.nanpercentile(arr, 95)),
        float(np.nanmax(arr)),
    )


def te_by_step(case_dir: Path) -> dict[int, tuple[float, float, float, float]]:
    out = {}
    for p in case_dir.glob("Te_out.*"):
        m = re.search(r"Te_out\.(\d+)$", p.name)
        if m:
            out[int(m.group(1))] = read_te_file(p)
    return out


def nearest_te(te: dict[int, tuple[float, float, float, float]], step: int) -> tuple[float, float, float, float]:
    if not te:
        return math.nan, math.nan, math.nan, math.nan
    key = step if step in te else min(te, key=lambda k: abs(k - step))
    return te[key]


def atom_temperature_from_ke(ke_eV: np.ndarray) -> np.ndarray:
    temp = np.full_like(ke_eV, np.nan, dtype=float)
    mobile = ke_eV > 1.0e-14
    temp[mobile] = 2.0 * ke_eV[mobile] / (3.0 * KB_EV_K)
    return temp


def neighbor_density_and_components(pos: np.ndarray, cutoff_A: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tree = cKDTree(pos)
    pairs = np.asarray(list(tree.query_pairs(cutoff_A)), dtype=int)
    n = len(pos)
    neigh = np.zeros(n, dtype=np.int32)
    parent = np.arange(n, dtype=np.int32)
    size = np.ones(n, dtype=np.int32)

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]

    if pairs.size:
        for a, b in pairs:
            neigh[a] += 1
            neigh[b] += 1
            union(int(a), int(b))

    roots = np.fromiter((find(i) for i in range(n)), dtype=np.int32, count=n)
    _, labels = np.unique(roots, return_inverse=True)
    comp_sizes = np.bincount(labels)
    return neigh.astype(float), labels.astype(np.int32), comp_sizes.astype(np.int32)


def nan_smooth(field: np.ndarray, sigma: float) -> np.ndarray:
    valid = np.isfinite(field)
    if not np.any(valid):
        return field
    vals = np.where(valid, field, 0.0)
    weights = valid.astype(float)
    sm_vals = gaussian_filter(vals, sigma=sigma, mode="nearest")
    sm_weights = gaussian_filter(weights, sigma=sigma, mode="nearest")
    out = np.full_like(field, np.nan, dtype=float)
    ok = sm_weights > 1.0e-8
    out[ok] = sm_vals[ok] / sm_weights[ok]
    return out


def binned_surface_p95(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    bounds: list[tuple[float, float]],
    mask: np.ndarray,
    nx: int = 8,
    ny: int = 8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xb = np.linspace(bounds[0][0], bounds[0][1], nx + 1)
    yb = np.linspace(bounds[1][0], bounds[1][1], ny + 1)
    surf = np.full((nx, ny), np.nan)
    for i in range(nx):
        xm = (x >= xb[i]) & (x < xb[i + 1])
        for j in range(ny):
            m = mask & xm & (y >= yb[j]) & (y < yb[j + 1])
            if np.count_nonzero(m) >= 4:
                surf[i, j] = float(np.nanpercentile(z[m], 95.0))
    return surf, 0.5 * (xb[:-1] + xb[1:]), 0.5 * (yb[:-1] + yb[1:])


def radial_bin_masks(xc: np.ndarray, yc: np.ndarray, bounds: list[tuple[float, float]]) -> tuple[np.ndarray, np.ndarray]:
    xmid = 0.5 * (bounds[0][0] + bounds[0][1])
    ymid = 0.5 * (bounds[1][0] + bounds[1][1])
    rr = np.sqrt((xc[:, None] - xmid) ** 2 + (yc[None, :] - ymid) ** 2)
    rmax = float(np.nanmax(rr))
    return rr <= 0.35 * rmax, rr >= 0.68 * rmax


def coarse_field(
    x_nm: np.ndarray,
    z_nm: np.ndarray,
    values: np.ndarray,
    x_edges: np.ndarray,
    z_edges: np.ndarray,
    reducer: str = "mean",
    min_count: int = 2,
) -> np.ndarray:
    out = np.full((len(z_edges) - 1, len(x_edges) - 1), np.nan)
    xi = np.digitize(x_nm, x_edges) - 1
    zi = np.digitize(z_nm, z_edges) - 1
    ok = (xi >= 0) & (xi < out.shape[1]) & (zi >= 0) & (zi < out.shape[0]) & np.isfinite(values)
    for iz in range(out.shape[0]):
        row = ok & (zi == iz)
        if not np.any(row):
            continue
        for ix in np.unique(xi[row]):
            m = row & (xi == ix)
            if np.count_nonzero(m) < min_count:
                continue
            if reducer == "p95":
                out[iz, ix] = float(np.nanpercentile(values[m], 95.0))
            elif reducer == "max":
                out[iz, ix] = float(np.nanmax(values[m]))
            else:
                out[iz, ix] = float(np.nanmean(values[m]))
    return out


def phase_code(temp: np.ndarray, vapor: np.ndarray, detached: np.ndarray, stable: np.ndarray, cluster: np.ndarray) -> np.ndarray:
    phase = np.zeros(len(temp), dtype=float)
    phase[temp > SI_MELT_K] = 1.0
    phase[vapor] = 2.0
    phase[detached] = 3.0
    phase[stable] = 4.0
    phase[cluster] = 5.0
    return phase


def fmt(x: float) -> str:
    if not math.isfinite(float(x)):
        return ""
    return f"{float(x):.8g}"


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[OK] {path}")


def analyze(project_root: Path, output_root: Path, case: str, equilibration_ps: float) -> None:
    case_dir = project_root / "results" / "v61_lammps_ttm_mod_ablation_scan" / case
    dump_path = case_dir / f"dump.v61_{case}.lammpstrj"
    log_path = case_dir / "log.lammps"
    if not dump_path.exists():
        raise FileNotFoundError(dump_path)

    out_light = output_root / "lightweight_results" / "v70_postprocess_diagnostics" / case
    out_fig = output_root / "figures" / "v70_postprocess_diagnostics" / case
    out_post = output_root / "post" / "v70_postprocess_diagnostics" / case
    out_light.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)
    out_post.mkdir(parents=True, exist_ok=True)

    thermo = parse_thermo(log_path)
    te = te_by_step(case_dir)
    target_times = np.asarray([5.0, 10.0, 20.0, 30.0, 50.0, 80.0, 100.0])
    best_fields: dict[float, dict[str, object]] = {}
    metrics: list[dict[str, str]] = []
    detached_history: dict[int, int] = {}
    prev_pos_by_id: dict[int, np.ndarray] | None = None
    initial = None
    initial_surface_bins = None
    initial_raw_surface_A = math.nan
    initial_surface_ref_A = math.nan
    rho0_neighbors = math.nan
    atom_volume_A3 = math.nan
    previous_time_ps = math.nan

    for iframe, (step, bounds, cols, data) in enumerate(iter_lammpstrj(dump_path)):
        idx = {name: i for i, name in enumerate(cols)}
        ids = data[:, idx["id"]].astype(int)
        x = data[:, idx["x"]]
        y = data[:, idx["y"]]
        z = data[:, idx["z"]]
        pos = np.column_stack([x, y, z])
        if all(k in idx for k in ("vx", "vy", "vz")):
            vx = data[:, idx["vx"]]
            vy = data[:, idx["vy"]]
            vz = data[:, idx["vz"]]
        elif prev_pos_by_id is not None and math.isfinite(previous_time_ps):
            current_time = thermo.get(step, {}).get("Time", step * 1.0e-4)
            dt = max(current_time - previous_time_ps, 1.0e-12)
            old = np.asarray([prev_pos_by_id.get(i, np.array([np.nan, np.nan, np.nan])) for i in ids])
            vel = (pos - old) / dt
            vx, vy, vz = vel[:, 0], vel[:, 1], vel[:, 2]
        else:
            vx = np.zeros_like(x)
            vy = np.zeros_like(x)
            vz = np.zeros_like(x)

        ke = data[:, idx["c_keatom"]] if "c_keatom" in idx else np.zeros_like(x)
        pe = data[:, idx["c_peatom"]] if "c_peatom" in idx else np.full_like(x, np.nan)
        temp = atom_temperature_from_ke(ke)
        if all(f"c_satom[{i}]" in idx for i in (1, 2, 3)):
            sxx = data[:, idx["c_satom[1]"]]
            syy = data[:, idx["c_satom[2]"]]
            szz = data[:, idx["c_satom[3]"]]
        else:
            sxx = syy = szz = np.zeros_like(x)

        neigh, comp_label, comp_sizes = neighbor_density_and_components(pos, cutoff_A=3.2)
        if iframe == 0:
            bulk_mask0 = (z > np.nanpercentile(z, 5)) & (z < np.nanpercentile(z, 85))
            rho0_neighbors = float(np.nanmean(neigh[bulk_mask0]))
            initial_raw_surface_A = float(np.nanmax(z))
            allmask = np.ones_like(z, dtype=bool)
            initial_surface_bins, xcent, ycent = binned_surface_p95(x, y, z, bounds, allmask)
            initial_surface_ref_A = float(np.nanmedian(initial_surface_bins))
            atom_volume_A3 = (
                (bounds[0][1] - bounds[0][0])
                * (bounds[1][1] - bounds[1][0])
                * max(initial_surface_ref_A - bounds[2][0], 1.0)
                / len(z)
            )
            initial = {
                "bounds": bounds,
                "surface_A": initial_surface_ref_A,
                "xcent": xcent,
                "ycent": ycent,
            }
        assert initial is not None
        assert initial_surface_bins is not None

        rho_ratio = neigh / max(rho0_neighbors, 1.0e-12)
        bulk_component = int(np.argmax(comp_sizes))
        in_bulk_component = comp_label == bulk_component
        surface_ref = initial["surface_A"]
        z_rel_nm = (z - surface_ref) * ANG_TO_NM
        xmid = 0.5 * (bounds[0][0] + bounds[0][1])
        ymid = 0.5 * (bounds[1][0] + bounds[1][1])
        x_center_nm = (x - xmid) * ANG_TO_NM
        y_center_nm = (y - ymid) * ANG_TO_NM

        vapor_conditions = np.column_stack([
            temp > SI_VAP_K,
            rho_ratio < 0.7,
            z > surface_ref + 5.0,
            vz > 0.0,
        ])
        vapor_candidate = np.sum(vapor_conditions, axis=1) >= 2
        detached_candidate = (z > surface_ref + 10.0) & (rho_ratio < 0.5) & (vz > 0.0)
        for atom_id, is_detached in zip(ids, detached_candidate):
            old = detached_history.get(int(atom_id), 0)
            detached_history[int(atom_id)] = old + 1 if is_detached else 0
        stable_by_history = np.asarray([detached_history.get(int(atom_id), 0) >= 3 for atom_id in ids])
        stable_ejecta = (z > surface_ref + 20.0) & (rho_ratio < 0.3) & (vz > 0.0) & stable_by_history

        detached_material = (~in_bulk_component) & (z > surface_ref)
        detached_components = [
            int(c)
            for c in np.unique(comp_label[detached_material])
            if np.count_nonzero(comp_label == c) > 0
        ]
        detached_cluster_sizes = [int(np.count_nonzero(comp_label == c)) for c in detached_components]
        cluster_labels = {c for c, s in zip(detached_components, detached_cluster_sizes) if s > 2}
        cluster_ejecta = np.asarray([label in cluster_labels for label in comp_label]) & (z > surface_ref)

        bad_for_surface = vapor_candidate | detached_candidate | stable_ejecta | (rho_ratio < 0.3)
        surface_mask = ~bad_for_surface
        surf_bins, xcent, ycent = binned_surface_p95(x, y, z, bounds, surface_mask)
        center_bins, edge_bins = radial_bin_masks(xcent, ycent, bounds)
        delta_bins = surf_bins - initial_surface_bins
        robust_lift_nm = float(np.nanmedian(delta_bins) * ANG_TO_NM)
        raw_lift_nm = float((np.nanmax(z) - initial_raw_surface_A) * ANG_TO_NM)
        surface_z_p95_global = float(np.nanpercentile(z[surface_mask], 95.0)) if np.any(surface_mask) else math.nan
        center_surface_A = float(np.nanmedian(surf_bins[center_bins])) if np.any(np.isfinite(surf_bins[center_bins])) else math.nan
        edge_surface_A = float(np.nanmedian(surf_bins[edge_bins])) if np.any(np.isfinite(surf_bins[edge_bins])) else math.nan
        center_initial_A = float(np.nanmedian(initial_surface_bins[center_bins]))
        edge_initial_A = float(np.nanmedian(initial_surface_bins[edge_bins]))
        crater_depth_center_nm = float(max(0.0, (edge_surface_A - center_surface_A) * ANG_TO_NM)) if math.isfinite(edge_surface_A) and math.isfinite(center_surface_A) else math.nan
        rim_lift_nm = float((edge_surface_A - edge_initial_A) * ANG_TO_NM) if math.isfinite(edge_surface_A) else math.nan

        near_surface = (z > surface_ref - 20.0) & (z < surface_ref + 15.0)
        rho_mean_surface = float(np.nanmean(rho_ratio[near_surface])) if np.any(near_surface) else math.nan
        rho_p10_surface = float(np.nanpercentile(rho_ratio[near_surface], 10.0)) if np.any(near_surface) else math.nan
        low_density_0p7 = int(np.count_nonzero(rho_ratio < 0.7))
        low_density_0p3 = int(np.count_nonzero(rho_ratio < 0.3))
        mean_vz_surface = float(np.nanmean(vz[near_surface])) if np.any(near_surface) else math.nan
        p95_vz_surface = float(np.nanpercentile(vz[near_surface], 95.0)) if np.any(near_surface) else math.nan
        upward_atoms = int(np.count_nonzero(near_surface & (vz > 0.0)))
        fast_upward_atoms = int(np.count_nonzero(near_surface & (vz > np.nanpercentile(vz[near_surface], 95.0)))) if np.any(near_surface) else 0
        outward_flux = float(np.nansum(vz[near_surface & (vz > 0.0)])) if np.any(near_surface) else math.nan

        pressure_proxy_GPa = -((sxx + syy + szz) / 3.0) / max(atom_volume_A3, 1.0e-12) * BAR_TO_GPA
        th = thermo.get(step, {})
        sim_time_ps = th.get("Time", step * 1.0e-4)
        time_ps = sim_time_ps - equilibration_ps
        te_min, te_mean, te_p95, te_max = nearest_te(te, step)
        surface_temp = temp[near_surface]
        temp_p95 = float(np.nanpercentile(surface_temp, 95.0)) if np.any(np.isfinite(surface_temp)) else math.nan
        temp_mean = float(np.nanmean(surface_temp)) if np.any(np.isfinite(surface_temp)) else math.nan

        row = {
            "case": case,
            "step": str(step),
            "time_ps": fmt(time_ps),
            "raw_lift_outlier_sensitive_nm": fmt(raw_lift_nm),
            "robust_lift_nm": fmt(robust_lift_nm),
            "surface_z_p95_global_nm": fmt((surface_z_p95_global - surface_ref) * ANG_TO_NM),
            "center_surface_z_p95_nm": fmt((center_surface_A - surface_ref) * ANG_TO_NM),
            "edge_surface_z_p95_nm": fmt((edge_surface_A - surface_ref) * ANG_TO_NM),
            "crater_depth_center_nm": fmt(crater_depth_center_nm),
            "rim_lift_nm": fmt(rim_lift_nm),
            "rho_mean_surface": fmt(rho_mean_surface),
            "rho_p10_surface": fmt(rho_p10_surface),
            "low_density_atoms_0p7": str(low_density_0p7),
            "low_density_atoms_0p3": str(low_density_0p3),
            "mean_vz_surface_A_ps": fmt(mean_vz_surface),
            "p95_vz_surface_A_ps": fmt(p95_vz_surface),
            "upward_atoms": str(upward_atoms),
            "fast_upward_atoms": str(fast_upward_atoms),
            "outward_flux_A_ps": fmt(outward_flux),
            "vapor_candidate_count": str(int(np.count_nonzero(vapor_candidate))),
            "detached_candidate_count": str(int(np.count_nonzero(detached_candidate))),
            "stable_ejecta_count": str(int(np.count_nonzero(stable_ejecta))),
            "cluster_ejecta_count": str(int(np.count_nonzero(cluster_ejecta))),
            "largest_ejecta_cluster_size": str(max(detached_cluster_sizes) if detached_cluster_sizes else 0),
            "ejected_mass_fraction": fmt(float(np.count_nonzero(stable_ejecta)) / len(ids)),
            "bulk_component_size": str(int(comp_sizes[bulk_component])),
            "detached_component_count": str(len(detached_components)),
            "detached_atom_count": str(int(np.count_nonzero(detached_material))),
            "largest_detached_cluster_size": str(max(detached_cluster_sizes) if detached_cluster_sizes else 0),
            "surface_atom_T_mean_K": fmt(temp_mean),
            "surface_atom_T_p95_K": fmt(temp_p95),
            "atom_T_max_proxy_K": fmt(float(np.nanmax(temp))),
            "melt_like_atom_count": str(int(np.count_nonzero(temp > SI_MELT_K))),
            "vapor_like_atom_count": str(int(np.count_nonzero(temp > SI_VAP_K))),
            "Te_min_K": fmt(te_min),
            "Te_mean_K": fmt(te_mean),
            "Te_p95_K": fmt(te_p95),
            "Te_max_K": fmt(te_max),
            "thermo_Tl_K": fmt(th.get("Temp", math.nan)),
            "thermo_pe_eV": fmt(th.get("PotEng", math.nan)),
            "thermo_ke_eV": fmt(th.get("KinEng", math.nan)),
            "thermo_total_eV": fmt(th.get("TotEng", math.nan)),
            "pressure_proxy_mean_GPa": fmt(float(np.nanmean(pressure_proxy_GPa))),
            "pressure_proxy_p01_GPa": fmt(float(np.nanpercentile(pressure_proxy_GPa, 1.0))),
            "pressure_proxy_p99_GPa": fmt(float(np.nanpercentile(pressure_proxy_GPa, 99.0))),
        }
        metrics.append(row)

        phase = phase_code(temp, vapor_candidate, detached_candidate, stable_ejecta, cluster_ejecta)
        for target in target_times:
            dist = abs(time_ps - target)
            old = best_fields.get(float(target))
            if old is None or dist < old["dist"]:
                y_section = np.abs(y_center_nm) <= 1.5
                x_edges = np.arange(-0.5 * (bounds[0][1] - bounds[0][0]) * ANG_TO_NM, 0.5 * (bounds[0][1] - bounds[0][0]) * ANG_TO_NM + 0.2, 0.2)
                z_edges = np.arange(-13.0, 6.2, 0.2)
                section = y_section
                count = coarse_field(x_center_nm[section], z_rel_nm[section], np.ones(np.count_nonzero(section)), x_edges, z_edges, reducer="mean", min_count=1)
                temp_field = coarse_field(x_center_nm[section], z_rel_nm[section], temp[section], x_edges, z_edges, reducer="mean")
                rho_field = coarse_field(x_center_nm[section], z_rel_nm[section], rho_ratio[section], x_edges, z_edges, reducer="mean")
                vz_field = coarse_field(x_center_nm[section], z_rel_nm[section], vz[section], x_edges, z_edges, reducer="mean")
                press_field = coarse_field(x_center_nm[section], z_rel_nm[section], pressure_proxy_GPa[section], x_edges, z_edges, reducer="mean")
                phase_field = coarse_field(x_center_nm[section], z_rel_nm[section], phase[section], x_edges, z_edges, reducer="max", min_count=1)
                best_fields[float(target)] = {
                    "dist": dist,
                    "target_time": float(target),
                    "time_ps": float(time_ps),
                    "step": int(step),
                    "x_edges": x_edges,
                    "z_edges": z_edges,
                    "count": count,
                    "temp": nan_smooth(temp_field, sigma=1.2),
                    "rho": nan_smooth(rho_field, sigma=1.2),
                    "vz": nan_smooth(vz_field, sigma=1.2),
                    "pressure": nan_smooth(press_field, sigma=1.2),
                    "phase": phase_field,
                    "scatter_x": x_center_nm[section],
                    "scatter_z": z_rel_nm[section],
                    "scatter_temp": np.clip(temp[section], 300.0, 12000.0),
                    "stable_ejecta": stable_ejecta[section],
                }

        prev_pos_by_id = {int(i): p.copy() for i, p in zip(ids, pos)}
        previous_time_ps = sim_time_ps

    write_csv(out_light / "ablation_metrics_physical.csv", metrics)
    write_field_npz(out_post / "stress_pressure_fields.npz", best_fields)
    plot_continuous_fields(out_fig / "si_ablation_continuous_fields_early.png", case, best_fields)
    plot_atom_scatter(out_fig / "atom_scatter_temperature_diagnostic.png", case, best_fields)
    plot_overview(out_fig / "overview_physical_diagnostics.png", case, metrics)
    write_readme(out_light / "README_V7_diagnostics.md", case, metrics)


def write_field_npz(path: Path, fields: dict[float, dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    for target, item in fields.items():
        tag = f"t{int(target):03d}ps"
        for key in ("temp", "rho", "vz", "pressure", "phase", "count", "x_edges", "z_edges"):
            arrays[f"{tag}_{key}"] = np.asarray(item[key])
        arrays[f"{tag}_actual_time_ps"] = np.asarray([item["time_ps"]], dtype=float)
        arrays[f"{tag}_step"] = np.asarray([item["step"]], dtype=int)
    np.savez_compressed(path, **arrays)
    print(f"[OK] {path}")


def plot_continuous_fields(path: Path, case: str, fields: dict[float, dict[str, object]]) -> None:
    ordered = [fields[t] for t in sorted(fields)]
    fig, axs = plt.subplots(4, len(ordered), figsize=(22, 10), sharex=True, sharey=True, constrained_layout=True)
    fig.suptitle(f"V7 continuous-field diagnostics for Si TTM-MD ablation: {case}", fontsize=18)
    row_defs = [
        ("temp", "atom T proxy / K", "turbo", 300.0, 8000.0),
        ("pressure", "pressure_proxy / GPa", "coolwarm", -0.8, 0.8),
        ("rho", "rho/rho0", "viridis", 0.0, 1.6),
        ("phase", "phase/ejecta code", "tab10", 0.0, 5.0),
    ]
    meshes = []
    for r, (key, label, cmap, vmin, vmax) in enumerate(row_defs):
        for c, item in enumerate(ordered):
            ax = axs[r, c]
            x_edges = np.asarray(item["x_edges"])
            z_edges = np.asarray(item["z_edges"])
            field = np.asarray(item[key])
            mesh = ax.pcolormesh(x_edges, z_edges, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            if c == len(ordered) - 1:
                meshes.append((mesh, label, r))
            ax.axhline(0.0, color="white", lw=0.7)
            ax.axhline(0.0, color="black", lw=0.35, ls="--")
            if r == 0:
                ax.set_title(f"{item['target_time']:.0f} ps\n(actual {item['time_ps']:.1f} ps)", fontsize=10)
            if c == 0:
                ax.set_ylabel("z from initial surface / nm")
            if r == 3:
                ax.set_xlabel("center x / nm")
            ax.set_xlim(float(x_edges[0]), float(x_edges[-1]))
            ax.set_ylim(float(z_edges[0]), float(z_edges[-1]))
    for mesh, label, r in meshes:
        cbar = fig.colorbar(mesh, ax=axs[r, :], shrink=0.78, pad=0.01)
        cbar.set_label(label)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def plot_atom_scatter(path: Path, case: str, fields: dict[float, dict[str, object]]) -> None:
    ordered = [fields[t] for t in sorted(fields)]
    fig, axes = plt.subplots(1, len(ordered), figsize=(22, 4.8), sharey=True, constrained_layout=True)
    fig.suptitle(f"V7 atom scatter temperature diagnostic: {case}", fontsize=16)
    sc = None
    for ax, item in zip(axes, ordered):
        sc = ax.scatter(
            np.asarray(item["scatter_x"]),
            np.asarray(item["scatter_z"]),
            c=np.asarray(item["scatter_temp"]),
            s=8,
            cmap="turbo",
            vmin=300,
            vmax=8000,
            linewidths=0,
            alpha=0.85,
        )
        stable = np.asarray(item["stable_ejecta"], dtype=bool)
        if np.any(stable):
            ax.scatter(np.asarray(item["scatter_x"])[stable], np.asarray(item["scatter_z"])[stable], s=28, facecolors="none", edgecolors="black", lw=0.8)
        ax.axhline(0.0, color="white", lw=0.8)
        ax.axhline(0.0, color="black", lw=0.4, ls="--")
        ax.set_title(f"{item['target_time']:.0f} ps", fontsize=10)
        ax.set_xlabel("center x / nm")
        ax.set_ylim(-13.0, 6.0)
    axes[0].set_ylabel("z from initial surface / nm")
    if sc is not None:
        cbar = fig.colorbar(sc, ax=axes, shrink=0.82, pad=0.01)
        cbar.set_label("atom temperature proxy / K")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def metric_series(rows: list[dict[str, str]], key: str) -> np.ndarray:
    out = []
    for row in rows:
        txt = row.get(key, "")
        out.append(float(txt) if txt else np.nan)
    return np.asarray(out, dtype=float)


def plot_overview(path: Path, case: str, rows: list[dict[str, str]]) -> None:
    t = metric_series(rows, "time_ps")
    fig, axs = plt.subplots(3, 3, figsize=(17, 12), constrained_layout=True)
    fig.suptitle(f"V7 physical diagnostics overview: {case}", fontsize=17)

    axs[0, 0].plot(t, metric_series(rows, "Te_p95_K"), label="Te p95", color="#d62728")
    axs[0, 0].plot(t, metric_series(rows, "thermo_Tl_K"), label="global Tl", color="#1f77b4")
    axs[0, 0].plot(t, metric_series(rows, "surface_atom_T_p95_K"), label="surface atom T p95", color="#ff7f0e")
    axs[0, 0].axhline(SI_MELT_K, color="k", ls="--", lw=0.9, label="Si melt")
    axs[0, 0].axhline(SI_VAP_K, color="0.4", ls=":", lw=0.9, label="Si vapor ref")
    axs[0, 0].set_ylabel("temperature / K")
    axs[0, 0].legend(fontsize=8)

    axs[0, 1].plot(t, metric_series(rows, "melt_like_atom_count"), label="Tproxy > melt")
    axs[0, 1].plot(t, metric_series(rows, "vapor_like_atom_count"), label="Tproxy > vapor")
    axs[0, 1].plot(t, metric_series(rows, "vapor_candidate_count"), label="vapor candidates")
    axs[0, 1].set_ylabel("atom count")
    axs[0, 1].legend(fontsize=8)

    axs[0, 2].plot(t, metric_series(rows, "raw_lift_outlier_sensitive_nm"), label="raw max-z lift", color="#d62728")
    axs[0, 2].plot(t, metric_series(rows, "robust_lift_nm"), label="robust p95 lift", color="#2ca02c")
    axs[0, 2].plot(t, metric_series(rows, "crater_depth_center_nm"), label="center crater depth", color="#9467bd")
    axs[0, 2].set_ylabel("nm")
    axs[0, 2].legend(fontsize=8)

    axs[1, 0].plot(t, metric_series(rows, "rho_mean_surface"), label="surface mean")
    axs[1, 0].plot(t, metric_series(rows, "rho_p10_surface"), label="surface p10")
    axs[1, 0].axhline(0.7, color="0.4", ls="--", lw=0.8)
    axs[1, 0].axhline(0.3, color="0.4", ls=":", lw=0.8)
    axs[1, 0].set_ylabel("rho/rho0")
    axs[1, 0].legend(fontsize=8)

    axs[1, 1].plot(t, metric_series(rows, "mean_vz_surface_A_ps"), label="mean vz surface")
    axs[1, 1].plot(t, metric_series(rows, "p95_vz_surface_A_ps"), label="p95 vz surface")
    axs[1, 1].set_ylabel("A/ps")
    axs[1, 1].legend(fontsize=8)

    axs[1, 2].plot(t, metric_series(rows, "detached_candidate_count"), label="detached candidate")
    axs[1, 2].plot(t, metric_series(rows, "stable_ejecta_count"), label="stable ejecta")
    axs[1, 2].plot(t, metric_series(rows, "cluster_ejecta_count"), label="cluster ejecta")
    axs[1, 2].set_ylabel("atom count")
    axs[1, 2].legend(fontsize=8)

    axs[2, 0].plot(t, metric_series(rows, "bulk_component_size"), label="bulk component")
    axs[2, 0].plot(t, metric_series(rows, "detached_atom_count"), label="detached atoms")
    axs[2, 0].set_ylabel("atom count")
    axs[2, 0].legend(fontsize=8)

    axs[2, 1].plot(t, metric_series(rows, "pressure_proxy_p99_GPa"), label="p99")
    axs[2, 1].plot(t, metric_series(rows, "pressure_proxy_p01_GPa"), label="p01")
    axs[2, 1].axhline(0.0, color="k", lw=0.8)
    axs[2, 1].set_ylabel("pressure_proxy / GPa")
    axs[2, 1].legend(fontsize=8)

    e0 = metric_series(rows, "thermo_total_eV")
    if np.isfinite(e0).any():
        axs[2, 2].plot(t, e0 - e0[np.where(np.isfinite(e0))[0][0]], label="Delta total")
    axs[2, 2].plot(t, metric_series(rows, "thermo_ke_eV"), label="atom KE")
    axs[2, 2].plot(t, metric_series(rows, "thermo_pe_eV"), label="atom PE")
    axs[2, 2].set_ylabel("eV")
    axs[2, 2].legend(fontsize=8)

    for ax in axs.ravel():
        ax.set_xlabel("time after laser start / ps")
        ax.grid(True, alpha=0.25)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def write_readme(path: Path, case: str, rows: list[dict[str, str]]) -> None:
    last = rows[-1]
    max_stable = int(np.nanmax(metric_series(rows, "stable_ejecta_count")))
    max_detached = int(np.nanmax(metric_series(rows, "detached_candidate_count")))
    max_vapor_cand = int(np.nanmax(metric_series(rows, "vapor_candidate_count")))
    max_vapor_like = int(np.nanmax(metric_series(rows, "vapor_like_atom_count")))
    max_raw = float(np.nanmax(metric_series(rows, "raw_lift_outlier_sensitive_nm")))
    max_robust = float(np.nanmax(metric_series(rows, "robust_lift_nm")))
    max_t = float(np.nanmax(metric_series(rows, "surface_atom_T_p95_K")))
    max_rho_low = int(np.nanmax(metric_series(rows, "low_density_atoms_0p3")))
    max_press = float(np.nanmax(metric_series(rows, "pressure_proxy_p99_GPa")))

    if max_stable <= 3 and max_detached < 20:
        state = "early overheating / melting / expansion, with only weak initial ejecta"
    elif max_stable < 100:
        state = "incipient ablation with limited stable ejecta"
    else:
        state = "clear atomistic ejection"

    lift_pollution = "yes" if max_raw > max_robust * 1.8 + 0.5 else "limited"
    text = f"""# V7 diagnostics: {case}

This file is generated from the existing V61 LAMMPS trajectory.  No simulation
parameters were changed.

## Main diagnosis

- Current state: **{state}**.
- Maximum surface atom temperature proxy p95: `{max_t:.3g} K`.
- Maximum vapor-like atom count by temperature proxy alone: `{max_vapor_like}`.
- Maximum vapor-candidate count after density/height/velocity filters: `{max_vapor_cand}`.
- Maximum detached-candidate count: `{max_detached}`.
- Maximum stable-ejecta count: `{max_stable}`.
- Maximum robust p95 surface lift: `{max_robust:.3g} nm`.
- Maximum raw max-z lift: `{max_raw:.3g} nm`.
- Raw lift outlier pollution: `{lift_pollution}`.
- Maximum low-density atoms with rho/rho0 < 0.3: `{max_rho_low}`.
- Maximum pressure proxy p99: `{max_press:.3g} GPa`.

## Why T > vapor can be large while stable ejecta is small

The per-atom temperature used here is a kinetic-energy proxy.  It can flag many
hot or strongly disordered atoms, but a real ejected atom must also be spatially
detached, low density, moving upward, and persistent across frames.  In this
trajectory, many atoms become hot, but most remain connected to the main slab.
That means the present 0.5 J/cm2, 100 ps run is not yet a strong plume case.

## Surface and crater metrics

`raw_lift_outlier_sensitive` uses the maximum z position and is kept only as a
diagnostic.  The main surface metric is `robust_lift`, computed from local
surface p95 values after excluding vapor/detached candidates and very low
density atoms.  Use `robust_lift`, `center_surface_z_p95`, `edge_surface_z_p95`,
`crater_depth_center`, and `rim_lift` for physical discussion.

## Pressure and stress

The pressure field is labelled `pressure_proxy`.  It is derived from per-atom
stress divided by an approximate initial atomic volume.  It is useful for
tracking stress-wave timing and sign, but should not be presented as a fully
calibrated hydrostatic pressure without additional validation.

## Is this enough for a paper-style main figure?

Not yet.  The continuous fields are useful diagnostics, but the current model is
too small laterally, only runs to 100 ps, and produces little stable ejecta.  V8
should update the physical model before treating the output as a main result.

## Likely V8 changes needed

- Increase total time to at least 300 ps.
- Use a larger lateral cell for pressure waves and plume formation.
- Re-check fluence normalization, reflected vs absorbed fluence, and optical
  penetration depth.
- Keep Si-specific TTM parameters explicit in the config.
- Add a controlled fluence scan rather than blindly increasing fluence.
- Treat Stillinger-Weber Si as a useful baseline but warn that high-temperature
  vapor/plume behavior requires validation.
"""
    path.write_text(text, encoding="utf-8")
    print(f"[OK] {path}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=Path.cwd())
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="fluence_0p5_100ps")
    ap.add_argument("--equilibration-ps", type=float, default=2.0)
    args = ap.parse_args()
    project_root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else project_root
    analyze(project_root, output_root, args.case, args.equilibration_ps)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
