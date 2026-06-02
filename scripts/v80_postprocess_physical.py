"""
v80_postprocess_physical.py
===========================

V8 corrected diagnostics for Si TTM-MD ablation trajectories.

The V7 diagnostics intentionally exposed problems in the earlier proxy fields.
This V8 postprocessor keeps the same input trajectory but corrects the main
diagnostic definitions:

* density field: x-z bin density normalized by the initial x-z bin density
* atom temperature: local-drift-subtracted thermal kinetic temperature
* stress/pressure: bin-summed LAMMPS stress/atom divided by bin volume
* surface: robust p95 profile excluding low-density and detached atoms

It is still a diagnostic for the existing trajectory; it does not rerun LAMMPS.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
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
AMU_A2_PS2_TO_EV = 1.0364269656262175e-4
BAR_TO_GPA = 1.0e-4
ANG_TO_NM = 0.1
SI_MELT_K = 1687.0
SI_VAP_K = 3538.0
HEATED_SURFACE_DEPTH_NM = 3.0
HEATED_SURFACE_ABOVE_NM = 0.5
MIN_FIELD_BIN_COUNT = 2
MIN_STRESS_BIN_COUNT = 3


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
                raise ValueError(bounds_header)
            bounds = []
            for _ in range(3):
                lo, hi, *_ = f.readline().split()
                bounds.append((float(lo), float(hi)))
            atom_header = f.readline().strip()
            if not atom_header.startswith("ITEM: ATOMS"):
                raise ValueError(atom_header)
            cols = atom_header.split()[2:]
            data = np.loadtxt([f.readline() for _ in range(natoms)])
            if data.ndim == 1:
                data = data.reshape(1, -1)
            yield step, bounds, cols, data


def iter_lammpstrj_many(paths: list[Path]):
    seen: set[int] = set()
    for path in paths:
        for step, bounds, cols, data in iter_lammpstrj(path):
            if step in seen:
                continue
            seen.add(step)
            yield step, bounds, cols, data


def parse_thermo(log_path: Path) -> dict[int, dict[str, float]]:
    out: dict[int, dict[str, float]] = {}
    if not log_path.exists():
        return out
    header = None
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
                out[int(row["Step"])] = row
    return out


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
    return float(np.nanmin(arr)), float(np.nanmean(arr)), float(np.nanpercentile(arr, 95)), float(np.nanmax(arr))


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


def smooth_nan(a: np.ndarray, sigma: float) -> np.ndarray:
    valid = np.isfinite(a)
    if not np.any(valid):
        return a
    values = np.where(valid, a, 0.0)
    weights = valid.astype(float)
    sv = gaussian_filter(values, sigma=sigma, mode="nearest")
    sw = gaussian_filter(weights, sigma=sigma, mode="nearest")
    out = np.full_like(a, np.nan, dtype=float)
    ok = sw > 1.0e-10
    out[ok] = sv[ok] / sw[ok]
    return out


def components_and_neighbor_density(pos: np.ndarray, cutoff_A: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tree = cKDTree(pos)
    pairs = np.asarray(list(tree.query_pairs(cutoff_A)), dtype=int)
    n = len(pos)
    counts = np.zeros(n, dtype=np.int32)
    parent = np.arange(n, dtype=np.int32)
    comp_size = np.ones(n, dtype=np.int32)

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if comp_size[ra] < comp_size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        comp_size[ra] += comp_size[rb]

    if pairs.size:
        for a, b in pairs:
            counts[a] += 1
            counts[b] += 1
            union(int(a), int(b))
    roots = np.fromiter((find(i) for i in range(n)), dtype=np.int32, count=n)
    _, labels = np.unique(roots, return_inverse=True)
    sizes = np.bincount(labels)
    return counts.astype(float), labels.astype(np.int32), sizes.astype(np.int32)


def edges_from_bounds(bounds: list[tuple[float, float]], initial_surface_A: float, dx_nm: float, dz_nm: float) -> tuple[np.ndarray, np.ndarray]:
    xw_nm = (bounds[0][1] - bounds[0][0]) * ANG_TO_NM
    x_edges = np.arange(-0.5 * xw_nm, 0.5 * xw_nm + dx_nm * 0.51, dx_nm)
    z_bottom_nm = (bounds[2][0] - initial_surface_A) * ANG_TO_NM
    z_edges = np.arange(math.floor(z_bottom_nm / dz_nm) * dz_nm, 8.0 + dz_nm * 0.51, dz_nm)
    return x_edges, z_edges


def bin_indices(x_nm: np.ndarray, z_nm: np.ndarray, x_edges: np.ndarray, z_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ix = np.digitize(x_nm, x_edges) - 1
    iz = np.digitize(z_nm, z_edges) - 1
    ok = (ix >= 0) & (ix < len(x_edges) - 1) & (iz >= 0) & (iz < len(z_edges) - 1)
    return ix, iz, ok


def mean_by_bin(values: np.ndarray, ix: np.ndarray, iz: np.ndarray, ok: np.ndarray, shape: tuple[int, int], min_count: int) -> np.ndarray:
    sumv = np.zeros(shape, dtype=float)
    count = np.zeros(shape, dtype=float)
    valid = ok & np.isfinite(values)
    np.add.at(sumv, (iz[valid], ix[valid]), values[valid])
    np.add.at(count, (iz[valid], ix[valid]), 1.0)
    out = np.full(shape, np.nan)
    m = count >= min_count
    out[m] = sumv[m] / count[m]
    return out


def max_by_bin(values: np.ndarray, ix: np.ndarray, iz: np.ndarray, ok: np.ndarray, shape: tuple[int, int], min_count: int) -> np.ndarray:
    out = np.full(shape, np.nan)
    count = np.zeros(shape, dtype=int)
    valid = ok & np.isfinite(values)
    for i, j, v in zip(ix[valid], iz[valid], values[valid]):
        if not np.isfinite(out[j, i]) or v > out[j, i]:
            out[j, i] = v
        count[j, i] += 1
    out[count < min_count] = np.nan
    return out


def count_by_bin(ix: np.ndarray, iz: np.ndarray, ok: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    out = np.zeros(shape, dtype=float)
    np.add.at(out, (iz[ok], ix[ok]), 1.0)
    return out


def local_thermal_temperature(
    mass_amu: float,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    ix: np.ndarray,
    iz: np.ndarray,
    ok: np.ndarray,
    shape: tuple[int, int],
) -> np.ndarray:
    vxm = mean_by_bin(vx, ix, iz, ok, shape, min_count=1)
    vym = mean_by_bin(vy, ix, iz, ok, shape, min_count=1)
    vzm = mean_by_bin(vz, ix, iz, ok, shape, min_count=1)
    lvx = np.zeros_like(vx)
    lvy = np.zeros_like(vy)
    lvz = np.zeros_like(vz)
    valid = ok.copy()
    if np.any(valid):
        valid_idx = np.where(valid)[0]
        finite_local = np.isfinite(vxm[iz[valid_idx], ix[valid_idx]])
        valid[valid_idx] = finite_local
    lvx[valid] = vxm[iz[valid], ix[valid]]
    lvy[valid] = vym[iz[valid], ix[valid]]
    lvz[valid] = vzm[iz[valid], ix[valid]]
    v2 = (vx - lvx) ** 2 + (vy - lvy) ** 2 + (vz - lvz) ** 2
    ke_th = 0.5 * mass_amu * v2 * AMU_A2_PS2_TO_EV
    temp = 2.0 * ke_th / (3.0 * KB_EV_K)
    temp[~valid] = np.nan
    return temp


def raw_temperature_from_ke(ke_eV: np.ndarray) -> np.ndarray:
    out = np.full_like(ke_eV, np.nan, dtype=float)
    mask = ke_eV > 1.0e-14
    out[mask] = 2.0 * ke_eV[mask] / (3.0 * KB_EV_K)
    return out


def virial_fields_gpa(
    sxx: np.ndarray,
    syy: np.ndarray,
    szz: np.ndarray,
    sxy: np.ndarray,
    sxz: np.ndarray,
    syz: np.ndarray,
    ix: np.ndarray,
    iz: np.ndarray,
    ok: np.ndarray,
    shape: tuple[int, int],
    bin_volume_A3: float,
    min_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    count = count_by_bin(ix, iz, ok, shape)
    fields = []
    for comp in (sxx, syy, szz, sxy, sxz, syz):
        acc = np.zeros(shape, dtype=float)
        valid = ok & np.isfinite(comp)
        np.add.at(acc, (iz[valid], ix[valid]), comp[valid])
        sigma = -acc / max(bin_volume_A3, 1.0e-30) * BAR_TO_GPA
        sigma[count < min_count] = np.nan
        fields.append(sigma)
    bx, by, bz, bxy, bxz, byz = fields
    pressure = (bx + by + bz) / 3.0
    vm = np.sqrt(0.5 * ((bx - by) ** 2 + (by - bz) ** 2 + (bz - bx) ** 2 + 6.0 * (bxy**2 + bxz**2 + byz**2)))
    return pressure, vm, smooth_nan(pressure, 1.0), smooth_nan(vm, 1.0)


def surface_profile(
    x_nm: np.ndarray,
    z_nm: np.ndarray,
    keep: np.ndarray,
    x_edges: np.ndarray,
    p: float = 95.0,
    min_count: int = 3,
) -> np.ndarray:
    out = np.full(len(x_edges) - 1, np.nan)
    ix = np.digitize(x_nm, x_edges) - 1
    for i in range(len(out)):
        m = keep & (ix == i)
        if np.count_nonzero(m) >= min_count:
            out[i] = float(np.nanpercentile(z_nm[m], p))
    return out


def fmt(x: float) -> str:
    if not math.isfinite(float(x)):
        return ""
    return f"{float(x):.8g}"


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[OK] {path}")


def resolve_case_inputs(project_root: Path, case: str) -> tuple[Path, list[Path], float]:
    v61_dir = project_root / "results" / "v61_lammps_ttm_mod_ablation_scan" / case
    v61_dump = v61_dir / f"dump.v61_{case}.lammpstrj"
    if v61_dump.exists():
        return v61_dir, [v61_dump], 2.0

    v80_dir = project_root / "results" / "v80_lammps_ttm_md_physical_model" / case
    v80_dumps = [
        v80_dir / f"dump.v80_{case}.early.lammpstrj",
        v80_dir / f"dump.v80_{case}.late.lammpstrj",
        v80_dir / f"dump.v80_{case}.lammpstrj",
    ]
    v80_dumps = [p for p in v80_dumps if p.exists()]
    if v80_dumps:
        return v80_dir, v80_dumps, 0.0

    raise FileNotFoundError(f"No V61/V80 dump found for case {case}")


def read_label_value_file(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    if not path.exists():
        return out
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for i in range(0, len(lines) - 1, 2):
        label = lines[i].strip()
        try:
            out[label] = float(lines[i + 1].strip())
        except ValueError:
            continue
    return out


def ttm_source_notes(case_dir: Path, case: str) -> str:
    candidates = list(case_dir.glob(f"v80_{case}.ttm_mod")) + list(case_dir.glob("*.ttm_mod"))
    if not candidates:
        return (
            "No `ttm_mod` parameter file was found for this case.  For a "
            "free-surface case, this is expected because `fix ttm/mod` is not used."
        )
    vals = read_label_value_file(candidates[0])
    i0 = vals.get("I_0, energy/(time*length^2) units", math.nan)
    lskin = vals.get("l_skin, length units", math.nan)
    tau = vals.get("tau, time units", math.nan)
    de = vals.get("D_e, length^2/time units", math.nan)
    return (
        f"`{candidates[0].name}` reports `I_0 = {i0:.6g}`, "
        f"`l_skin = {lskin:.6g} A` ({lskin * ANG_TO_NM:.3g} nm), "
        f"`tau = {tau:.6g} ps`, and `D_e = {de:.6g}`.  The source is not "
        "intended as uniform slab heating, but the current penetration depth and "
        "electron diffusion parameters can still make the 30 ps diagnostic look "
        "too deep or too spatially uniform.  Treat this as a parameter-validation "
        "warning before any paper-style Si interpretation."
    )


def case_scope_label(case: str) -> str:
    if "ttm_periodic_slab_diagnostic" in case:
        return "TTM-MD periodic slab diagnostic, 30 ps; not strict free-surface ablation"
    if "free_surface_ablation" in case:
        return "Free-surface ablation morphology case; non-periodic z, no fix ttm/mod"
    return "V8 physical diagnostics"


def analyze(project_root: Path, output_root: Path, case: str, equilibration_ps: float, mass_amu: float) -> None:
    case_dir, dump_paths, inferred_time_offset_ps = resolve_case_inputs(project_root, case)
    if equilibration_ps < 0:
        equilibration_ps = inferred_time_offset_ps
    log_path = case_dir / "log.lammps"

    out_light = output_root / "lightweight_results" / "v80_physical_model_update" / case
    out_fig = output_root / "figures" / "v80_physical_model_update" / case
    out_post = output_root / "post" / "v80_physical_model_update" / case
    out_light.mkdir(parents=True, exist_ok=True)
    out_fig.mkdir(parents=True, exist_ok=True)
    out_post.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] input dir : {case_dir}")
    print(f"[INFO] figure dir: {out_fig}")
    print(f"[INFO] scope     : {case_scope_label(case)}")

    thermo = parse_thermo(log_path)
    te = te_by_step(case_dir)
    if "20ps_free_surface_ablation_debug" in case:
        target_early = [0, 10, 12, 14, 15, 16, 18, 20]
        target_ejecta = [0, 10, 12, 14, 15, 16, 18, 20]
    elif "50ps_free_surface_ablation_stable_debug" in case:
        target_early = [0, 10, 20, 30, 40, 50]
        target_ejecta = [0, 10, 20, 30, 40, 50]
    elif "free_surface_ablation_quickcheck" in case:
        target_early = [0, 10, 30, 60, 100]
        target_ejecta = [0, 10, 30, 60, 100]
    elif "free_surface_ablation" in case:
        target_early = [0, 10, 30, 100, 300]
        target_ejecta = [0, 10, 30, 100, 300]
    elif "ttm_periodic_slab_diagnostic" in case:
        target_early = [0, 5, 10, 20, 30]
        target_ejecta = [0, 5, 10, 20, 30]
    else:
        target_early = [5, 10, 20, 30, 50, 80, 100]
        target_ejecta = [30, 60, 100]
    targets = sorted(set(target_early + target_ejecta))
    metrics: list[dict[str, str]] = []
    surface_rows: list[dict[str, str]] = []
    heat_group_rows: list[dict[str, str]] = []
    fields: dict[int, dict[str, object]] = {}
    detached_history: dict[int, int] = {}

    initial_surface_A = None
    initial_count_field = None
    initial_surface_profile = None
    x_edges = z_edges = None
    x_centers = None
    rho0_neighbor = math.nan
    bin_volume_A3 = math.nan

    for iframe, (step, bounds, cols, data) in enumerate(iter_lammpstrj_many(dump_paths)):
        idx = {name: i for i, name in enumerate(cols)}
        ids = data[:, idx["id"]].astype(int)
        x = data[:, idx["x"]]
        y = data[:, idx["y"]]
        z = data[:, idx["z"]]
        vx = data[:, idx["vx"]]
        vy = data[:, idx["vy"]]
        vz = data[:, idx["vz"]]
        ke = data[:, idx["c_keatom"]]
        sxx = data[:, idx["c_satom[1]"]]
        syy = data[:, idx["c_satom[2]"]]
        szz = data[:, idx["c_satom[3]"]]
        sxy = data[:, idx["c_satom[4]"]]
        sxz = data[:, idx["c_satom[5]"]]
        syz = data[:, idx["c_satom[6]"]]

        th = thermo.get(step, {})
        sim_time_ps = th.get("Time", step * 1.0e-4)
        time_ps = sim_time_ps - equilibration_ps

        if initial_surface_A is None:
            initial_surface_A = float(np.nanpercentile(z, 99.5))
            x_edges, z_edges = edges_from_bounds(bounds, initial_surface_A, dx_nm=0.2, dz_nm=0.2)
            x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])

        assert initial_surface_A is not None
        assert x_edges is not None and z_edges is not None and x_centers is not None

        xmid = 0.5 * (bounds[0][0] + bounds[0][1])
        ymid = 0.5 * (bounds[1][0] + bounds[1][1])
        x_nm = (x - xmid) * ANG_TO_NM
        z_nm = (z - initial_surface_A) * ANG_TO_NM
        y_nm = (y - ymid) * ANG_TO_NM
        # V8 continuous fields use a full-y x-z projection.  The earlier thin
        # center slice made pressure and bin-density fields too sparse for a
        # small debug cell.
        section = np.ones_like(y_nm, dtype=bool)
        shape = (len(z_edges) - 1, len(x_edges) - 1)
        ix, iz, ok_all = bin_indices(x_nm, z_nm, x_edges, z_edges)
        ok_section = ok_all & section
        if (not math.isfinite(bin_volume_A3)) or bin_volume_A3 <= 0:
            dx_A = (x_edges[1] - x_edges[0]) / ANG_TO_NM
            dz_A = (z_edges[1] - z_edges[0]) / ANG_TO_NM
            dy_A = bounds[1][1] - bounds[1][0]
            bin_volume_A3 = dx_A * dy_A * dz_A

        raw_temp = raw_temperature_from_ke(ke)
        thermal_temp = local_thermal_temperature(mass_amu, vx, vy, vz, ix, iz, ok_all, shape)
        neigh, comp_label, comp_sizes = components_and_neighbor_density(np.column_stack([x, y, z]), cutoff_A=3.2)
        if iframe == 0:
            bulk0 = (z_nm > -10.0) & (z_nm < -1.0)
            rho0_neighbor = float(np.nanmean(neigh[bulk0]))
            initial_count_field = count_by_bin(ix[section], iz[section], ok_section[section], shape)
            initial_surface_profile = surface_profile(x_nm, z_nm, np.ones_like(z, dtype=bool), x_edges)
            if not np.nanmean(initial_count_field[initial_count_field > 0]) > 0:
                raise RuntimeError("Initial density field is empty")

        assert initial_count_field is not None
        assert initial_surface_profile is not None

        bin_count = count_by_bin(ix[section], iz[section], ok_section[section], shape)
        density_ratio = np.full(shape, np.nan)
        valid_initial = initial_count_field >= 2
        density_ratio[valid_initial] = bin_count[valid_initial] / initial_count_field[valid_initial]
        low_density_area_fraction = float(np.count_nonzero((density_ratio < 0.7) & valid_initial) / np.count_nonzero(valid_initial))
        void_area_fraction = float(np.count_nonzero((density_ratio < 0.3) & valid_initial) / np.count_nonzero(valid_initial))
        initial_density_mean_check = float(np.nanmean(density_ratio[valid_initial])) if iframe == 0 else math.nan

        neighbor_density = neigh / max(rho0_neighbor, 1.0e-12)
        bulk_label = int(np.argmax(comp_sizes))
        in_bulk = comp_label == bulk_label
        vapor_conditions = np.column_stack([
            thermal_temp > SI_VAP_K,
            neighbor_density < 0.7,
            z_nm > 0.5,
            vz > 0.0,
        ])
        vapor_candidate = np.sum(vapor_conditions, axis=1) >= 2
        detached_candidate = (z_nm > 1.0) & (neighbor_density < 0.5) & (vz > 0.0)
        for atom_id, is_detached in zip(ids, detached_candidate):
            detached_history[int(atom_id)] = detached_history.get(int(atom_id), 0) + 1 if is_detached else 0
        stable_history = np.asarray([detached_history.get(int(atom_id), 0) >= 3 for atom_id in ids])
        stable_ejecta = (z_nm > 2.0) & (vz > 0.0) & ((neighbor_density < 0.3) | (~in_bulk)) & stable_history

        detached_material = (~in_bulk) & (z_nm > 0.0)
        detached_labels = [int(label) for label in np.unique(comp_label[detached_material])]
        detached_sizes = [int(np.count_nonzero(comp_label == label)) for label in detached_labels]
        cluster_labels = {label for label, size in zip(detached_labels, detached_sizes) if size > 2}
        cluster_ejecta = np.asarray([label in cluster_labels for label in comp_label]) & (z_nm > 0.0)

        keep_surface = ~(vapor_candidate | detached_candidate | stable_ejecta | (neighbor_density < 0.3))
        surf = surface_profile(x_nm, z_nm, keep_surface, x_edges)
        surf_delta = surf - initial_surface_profile
        robust_lift = float(np.nanmedian(surf_delta))
        raw_lift = float(np.nanmax(z_nm))
        center = np.abs(x_centers) <= 0.25 * max(abs(x_edges[0]), abs(x_edges[-1]))
        edge = np.abs(x_centers) >= 0.65 * max(abs(x_edges[0]), abs(x_edges[-1]))
        center_surf = float(np.nanmedian(surf[center]))
        edge_surf = float(np.nanmedian(surf[edge]))
        crater_depth = float(max(0.0, edge_surf - center_surf)) if math.isfinite(edge_surf) and math.isfinite(center_surf) else math.nan
        rim_lift = float(np.nanmedian(surf_delta[edge])) if np.any(np.isfinite(surf_delta[edge])) else math.nan

        pressure, vm_stress, pressure_sm, vm_stress_sm = virial_fields_gpa(
            sxx[section],
            syy[section],
            szz[section],
            sxy[section],
            sxz[section],
            syz[section],
            ix[section],
            iz[section],
            ok_section[section],
            shape,
            bin_volume_A3,
            min_count=MIN_STRESS_BIN_COUNT,
        )
        temp_field = mean_by_bin(thermal_temp[section], ix[section], iz[section], ok_section[section], shape, min_count=MIN_FIELD_BIN_COUNT)
        raw_temp_field = mean_by_bin(raw_temp[section], ix[section], iz[section], ok_section[section], shape, min_count=MIN_FIELD_BIN_COUNT)
        vz_field = mean_by_bin(vz[section], ix[section], iz[section], ok_section[section], shape, min_count=MIN_FIELD_BIN_COUNT)
        phase = np.zeros_like(thermal_temp)
        phase[thermal_temp > SI_MELT_K] = 1
        phase[vapor_candidate] = 2
        phase[detached_candidate] = 3
        phase[stable_ejecta] = 4
        phase[cluster_ejecta] = 5
        phase_field = max_by_bin(phase[section], ix[section], iz[section], ok_section[section], shape, min_count=1)

        surface_region = (z_nm >= -HEATED_SURFACE_DEPTH_NM) & (z_nm <= HEATED_SURFACE_ABOVE_NM) & keep_surface
        ejected_region = stable_ejecta
        surface_thermal = thermal_temp[surface_region]
        surface_raw = raw_temp[surface_region]
        finite_thermal = thermal_temp[np.isfinite(thermal_temp)]
        finite_ke = ke[np.isfinite(ke)]
        max_ke_eV = float(np.nanmax(finite_ke)) if finite_ke.size else math.nan
        max_vz_A_ps = float(np.nanmax(vz)) if vz.size else math.nan
        zmax_nm = float(np.nanmax(z_nm)) if z_nm.size else math.nan
        if finite_ke.size:
            top_ke_cut = float(np.nanpercentile(finite_ke, 99))
            top_ke = ke >= top_ke_cut
            top1_zmin_nm = float(np.nanmin(z_nm[top_ke])) if np.any(top_ke) else math.nan
            top1_zmax_nm = float(np.nanmax(z_nm[top_ke])) if np.any(top_ke) else math.nan
        else:
            top_ke_cut = top1_zmin_nm = top1_zmax_nm = math.nan
        te_min, te_mean, te_p95, te_max = nearest_te(te, step)
        row = {
            "case": case,
            "step": str(step),
            "time_ps": fmt(time_ps),
            "max_atom_ke_eV": fmt(max_ke_eV),
            "max_vz_A_ps": fmt(max_vz_A_ps),
            "surface_zmax_nm": fmt(zmax_nm),
            "top1pct_ke_threshold_eV": fmt(top_ke_cut),
            "top1pct_ke_zmin_nm": fmt(top1_zmin_nm),
            "top1pct_ke_zmax_nm": fmt(top1_zmax_nm),
            "raw_lift_outlier_sensitive_nm": fmt(raw_lift),
            "robust_p95_lift_nm": fmt(robust_lift),
            "center_crater_depth_nm": fmt(crater_depth),
            "rim_lift_nm": fmt(rim_lift),
            "center_surface_z_p95_nm": fmt(center_surf),
            "edge_surface_z_p95_nm": fmt(edge_surf),
            "bin_density_mean_check_initial": fmt(initial_density_mean_check),
            "low_density_area_fraction": fmt(low_density_area_fraction),
            "void_area_fraction": fmt(void_area_fraction),
            "heated_surface_region_depth_nm": fmt(HEATED_SURFACE_DEPTH_NM),
            "heated_surface_region_above_surface_nm": fmt(HEATED_SURFACE_ABOVE_NM),
            "heated_surface_atom_count": str(int(np.count_nonzero(surface_region))),
            "continuous_field_min_atom_count_per_bin": str(MIN_FIELD_BIN_COUNT),
            "stress_field_min_atom_count_per_bin": str(MIN_STRESS_BIN_COUNT),
            "neighbor_density_surface_mean": fmt(float(np.nanmean(neighbor_density[surface_region])) if np.any(surface_region) else math.nan),
            "atom_T_raw_mean_K": fmt(float(np.nanmean(surface_raw)) if surface_raw.size else math.nan),
            "atom_T_raw_p95_K": fmt(float(np.nanpercentile(surface_raw, 95)) if surface_raw.size else math.nan),
            "atom_T_raw_max_K": fmt(float(np.nanmax(surface_raw)) if surface_raw.size else math.nan),
            "atom_T_thermal_mean_K": fmt(float(np.nanmean(surface_thermal)) if surface_thermal.size else math.nan),
            "atom_T_thermal_p95_K": fmt(float(np.nanpercentile(surface_thermal, 95)) if surface_thermal.size else math.nan),
            "atom_T_thermal_max_K": fmt(float(np.nanmax(surface_thermal)) if surface_thermal.size else math.nan),
            "global_atom_T_thermal_mean_K": fmt(float(np.nanmean(finite_thermal)) if finite_thermal.size else math.nan),
            "molten_atom_count": str(int(np.count_nonzero(thermal_temp > SI_MELT_K))),
            "molten_atom_fraction": fmt(float(np.count_nonzero(thermal_temp > SI_MELT_K)) / len(ids)),
            "vapor_threshold_atom_count": str(int(np.count_nonzero(thermal_temp > SI_VAP_K))),
            "vapor_threshold_atom_fraction": fmt(float(np.count_nonzero(thermal_temp > SI_VAP_K)) / len(ids)),
            "ejecta_T_thermal_mean_K": fmt(float(np.nanmean(thermal_temp[ejected_region])) if np.any(ejected_region) else math.nan),
            "Te_min_K": fmt(te_min),
            "Te_mean_K": fmt(te_mean),
            "Te_p95_K": fmt(te_p95),
            "Te_max_K": fmt(te_max),
            "thermo_Tl_K": fmt(th.get("Temp", math.nan)),
            "upward_atoms": str(int(np.count_nonzero(surface_region & (vz > 0.0)))),
            "vapor_candidate_count": str(int(np.count_nonzero(vapor_candidate))),
            "detached_candidate_count": str(int(np.count_nonzero(detached_candidate))),
            "stable_ejecta_count": str(int(np.count_nonzero(stable_ejecta))),
            "cluster_ejecta_count": str(int(np.count_nonzero(cluster_ejecta))),
            "largest_ejecta_cluster_size": str(max(detached_sizes) if detached_sizes else 0),
            "ejected_mass_fraction": fmt(float(np.count_nonzero(stable_ejecta)) / len(ids)),
            "bulk_component_size": str(int(comp_sizes[bulk_label])),
            "bulk_component_fraction": fmt(float(comp_sizes[bulk_label]) / len(ids)),
            "bulk_component_delta_count": str(int(comp_sizes[bulk_label]) - len(ids)),
            "detached_atom_count": str(int(np.count_nonzero(detached_material))),
            "pressure_GPa_p01": fmt(float(np.nanpercentile(pressure, 1)) if np.any(np.isfinite(pressure)) else math.nan),
            "pressure_GPa_p99": fmt(float(np.nanpercentile(pressure, 99)) if np.any(np.isfinite(pressure)) else math.nan),
            "von_mises_stress_GPa_p99": fmt(float(np.nanpercentile(vm_stress, 99)) if np.any(np.isfinite(vm_stress)) else math.nan),
            "thermo_pe_eV": fmt(th.get("PotEng", math.nan)),
            "thermo_ke_eV": fmt(th.get("KinEng", math.nan)),
            "thermo_total_eV": fmt(th.get("TotEng", math.nan)),
            "thermo_max_atom_ke_eV": fmt(th.get("c_max_ke", math.nan)),
            "thermo_max_vz_A_ps": fmt(th.get("c_max_vz", math.nan)),
            "thermo_zmax_A": fmt(th.get("c_zmax", math.nan)),
        }
        metrics.append(row)
        for xc, zs in zip(x_centers, surf):
            surface_rows.append({"case": case, "step": str(step), "time_ps": fmt(time_ps), "x_nm": fmt(float(xc)), "surface_z_nm": fmt(float(zs))})

        for layer_i in range(1, 11):
            lo_nm = float(layer_i - 1)
            hi_nm = float(layer_i)
            layer_mask = (z_nm <= -lo_nm) & (z_nm > -hi_nm)
            temp_layer = thermal_temp[layer_mask]
            ke_layer = ke[layer_mask]
            vz_layer = vz[layer_mask]
            heat_group_rows.append({
                "case": case,
                "step": str(step),
                "time_ps": fmt(time_ps),
                "heat_group": f"heat_{layer_i:02d}",
                "depth_lo_nm": fmt(lo_nm),
                "depth_hi_nm": fmt(hi_nm),
                "atom_count": str(int(np.count_nonzero(layer_mask))),
                "mean_T_thermal_K": fmt(float(np.nanmean(temp_layer)) if temp_layer.size else math.nan),
                "p95_T_thermal_K": fmt(float(np.nanpercentile(temp_layer, 95)) if temp_layer.size else math.nan),
                "max_ke_eV": fmt(float(np.nanmax(ke_layer)) if ke_layer.size else math.nan),
                "max_vz_A_ps": fmt(float(np.nanmax(vz_layer)) if vz_layer.size else math.nan),
                "zmin_nm": fmt(float(np.nanmin(z_nm[layer_mask])) if np.any(layer_mask) else math.nan),
                "zmax_nm": fmt(float(np.nanmax(z_nm[layer_mask])) if np.any(layer_mask) else math.nan),
            })

        for target in targets:
            old = fields.get(target)
            dist = abs(time_ps - target)
            if old is None or dist < old["dist"]:
                fields[target] = {
                    "dist": dist,
                    "target_time_ps": target,
                    "actual_time_ps": time_ps,
                    "step": step,
                    "x_edges": x_edges.copy(),
                    "z_edges": z_edges.copy(),
                    "temperature": smooth_nan(temp_field, 1.0),
                    "temperature_raw": smooth_nan(raw_temp_field, 1.0),
                    "density": smooth_nan(density_ratio, 1.0),
                    "pressure": pressure_sm,
                    "stress_vm": vm_stress_sm,
                    "phase": phase_field,
                    "vz": smooth_nan(vz_field, 1.0),
                    "surface_x": x_centers.copy(),
                    "surface_z": surf.copy(),
                }

    write_csv(out_light / "ablation_metrics_physical_v8.csv", metrics)
    write_csv(out_light / "surface_profiles_v8.csv", surface_rows)
    write_csv(out_light / "heat_group_diagnostics_v8.csv", heat_group_rows)
    write_npz(out_post / "stress_pressure_fields_v8.npz", fields)
    plot_field_grid(out_fig / "si_ablation_continuous_fields_early_v8.png", case, fields, target_early, mode="early")
    plot_field_grid(out_fig / "si_ablation_continuous_fields_ejecta_v8.png", case, fields, target_ejecta, mode="ejecta")
    plot_overview(out_fig / "overview_physical_v8.png", case, metrics)
    plot_surface_profiles(out_fig / "surface_profile_evolution.png", case, surface_rows, target_early)
    plot_temperature_time(out_fig / "temperature_lattice_vs_time_v8.png", case, metrics)
    plot_temperature_time(out_fig / "temperature_diagnostics_v8.png", case, metrics)
    plot_pressure_time(out_fig / "pressure_vs_time_v8.png", case, metrics)
    plot_stability_diagnostics(out_fig / "stability_diagnostics_v8.png", case, metrics)
    plot_depth_temperature_profiles(out_fig / "depth_resolved_temperature_profile_v8.png", case, fields, target_early)
    if np.any(np.isfinite(series(metrics, "Te_mean_K"))):
        plot_electron_temperature_time(out_fig / "electron_temperature_vs_time_v8.png", case, metrics)
    write_readme(out_light / "README_V8_physical_model.md", case, metrics, case_dir)


def write_npz(path: Path, fields: dict[int, dict[str, object]]) -> None:
    arrays = {}
    for target, item in fields.items():
        tag = f"t{target:03d}ps"
        for key in ("x_edges", "z_edges", "temperature", "temperature_raw", "density", "pressure", "stress_vm", "phase", "vz", "surface_x", "surface_z"):
            arrays[f"{tag}_{key}"] = np.asarray(item[key])
        arrays[f"{tag}_actual_time_ps"] = np.asarray([item["actual_time_ps"]])
        arrays[f"{tag}_step"] = np.asarray([item["step"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays)
    print(f"[OK] {path}")


def percentile_limits(fields: list[np.ndarray], lo: float, hi: float) -> tuple[float, float]:
    chunks = [f[np.isfinite(f)].ravel() for f in fields if np.any(np.isfinite(f))]
    if not chunks:
        return 0.0, 1.0
    vals = np.concatenate(chunks)
    if vals.size == 0:
        return 0.0, 1.0
    return float(np.nanpercentile(vals, lo)), float(np.nanpercentile(vals, hi))


def plot_field_grid(path: Path, case: str, fields: dict[int, dict[str, object]], targets: list[int], mode: str) -> None:
    ordered = [fields[t] for t in targets if t in fields]
    if mode == "early":
        rows = [
            ("temperature", "thermal atom T / K", "turbo", 300.0, 8000.0),
            ("pressure", "hydrostatic pressure / GPa", "coolwarm", None, None),
            ("density", "bin density rho/rho0", "viridis", 0.0, 1.6),
            ("phase", "phase/ejecta code", "tab10", 0.0, 5.0),
        ]
        title = f"{case_scope_label(case)}\ncorrected continuous fields, early stage: {case}"
    else:
        rows = [
            ("temperature", "thermal atom T / K", "turbo", 300.0, 8000.0),
            ("density", "bin density rho/rho0", "viridis", 0.0, 1.6),
            ("phase", "phase/ejecta code", "tab10", 0.0, 5.0),
            ("vz", "vz / A ps^-1", "coolwarm", None, None),
        ]
        title = f"{case_scope_label(case)}\ncorrected continuous fields, ejecta window: {case}"
    fig, axs = plt.subplots(len(rows), len(ordered), figsize=(3.2 * len(ordered) + 3.0, 10), sharex=True, sharey=True, constrained_layout=True)
    if len(ordered) == 1:
        axs = axs[:, None]
    fig.suptitle(title, fontsize=16)
    for r, (key, label, cmap, vmin, vmax) in enumerate(rows):
        if vmin is None:
            vmin, vmax = percentile_limits([np.asarray(item[key]) for item in ordered], 1, 99)
            if abs(vmin - vmax) < 1.0e-12:
                vmax = vmin + 1.0
        last_mesh = None
        for c, item in enumerate(ordered):
            ax = axs[r, c]
            mesh = ax.pcolormesh(item["x_edges"], item["z_edges"], item[key], shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            last_mesh = mesh
            ax.plot(item["surface_x"], item["surface_z"], color="black", lw=0.9)
            ax.axhline(0.0, color="white", lw=0.7)
            ax.axhline(0.0, color="black", lw=0.35, ls="--")
            if r == 0:
                ax.set_title(f"{item['target_time_ps']:.0f} ps\nactual {item['actual_time_ps']:.1f} ps", fontsize=10)
            if c == 0:
                ax.set_ylabel("z from initial surface / nm")
            if r == len(rows) - 1:
                ax.set_xlabel("center x / nm")
            ax.set_ylim(-13.0, 8.0)
        if last_mesh is not None:
            cbar = fig.colorbar(last_mesh, ax=axs[r, :], shrink=0.78, pad=0.01)
            cbar.set_label(label)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def series(rows: list[dict[str, str]], key: str) -> np.ndarray:
    out = []
    for row in rows:
        try:
            out.append(float(row.get(key, "")))
        except ValueError:
            out.append(np.nan)
    return np.asarray(out, dtype=float)


def plot_temperature_time(path: Path, case: str, rows: list[dict[str, str]]) -> None:
    t = series(rows, "time_ps")
    fig, ax = plt.subplots(figsize=(9, 5.4), constrained_layout=True)
    fig.suptitle(f"{case_scope_label(case)}\ntemperature diagnostics: {case}", fontsize=13)
    curves = [
        ("thermo_Tl_K", "global lattice/atom T from thermo"),
        ("atom_T_thermal_mean_K", "heated surface thermal T mean"),
        ("atom_T_thermal_p95_K", "heated surface thermal T p95"),
        ("atom_T_thermal_max_K", "heated surface thermal T max - hotspot reference only"),
        ("global_atom_T_thermal_mean_K", "global drift-removed atom T mean"),
    ]
    for key, label in curves:
        y = series(rows, key)
        if np.any(np.isfinite(y)):
            ax.plot(t, y, marker="o", ms=3, label=label)
    ax.axhline(SI_MELT_K, color="k", ls="--", lw=0.9, label="Si melt 1687 K")
    ax.axhline(SI_VAP_K, color="0.35", ls=":", lw=0.9, label="vapor reference 3538 K")
    ax.set_xlabel("time / ps")
    ax.set_ylabel("temperature / K")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def plot_electron_temperature_time(path: Path, case: str, rows: list[dict[str, str]]) -> None:
    t = series(rows, "time_ps")
    fig, ax = plt.subplots(figsize=(8.6, 5.0), constrained_layout=True)
    fig.suptitle(f"{case_scope_label(case)}\nelectron temperature diagnostics: {case}", fontsize=13)
    for key, label in (
        ("Te_mean_K", "Te mean"),
        ("Te_p95_K", "Te p95"),
        ("Te_max_K", "Te max"),
    ):
        y = series(rows, key)
        if np.any(np.isfinite(y)):
            ax.plot(t, y, marker="o", ms=3, label=label)
    ax.set_xlabel("time / ps")
    ax.set_ylabel("electron temperature / K")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def plot_pressure_time(path: Path, case: str, rows: list[dict[str, str]]) -> None:
    t = series(rows, "time_ps")
    fig, ax = plt.subplots(figsize=(8.6, 5.0), constrained_layout=True)
    fig.suptitle(f"{case_scope_label(case)}\npressure/stress diagnostics: {case}", fontsize=13)
    for key, label in (
        ("pressure_GPa_p99", "hydrostatic pressure p99"),
        ("pressure_GPa_p01", "hydrostatic pressure p01"),
        ("von_mises_stress_GPa_p99", "von Mises stress p99"),
    ):
        y = series(rows, key)
        if np.any(np.isfinite(y)):
            ax.plot(t, y, marker="o", ms=3, label=label)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel("time / ps")
    ax.set_ylabel("GPa")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def plot_stability_diagnostics(path: Path, case: str, rows: list[dict[str, str]]) -> None:
    t = series(rows, "time_ps")
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.2), constrained_layout=True)
    fig.suptitle(f"{case_scope_label(case)}\npre-lost-atoms stability diagnostics: {case}", fontsize=13)

    axs[0, 0].plot(t, series(rows, "max_atom_ke_eV"), marker="o", ms=3, label="dump max atom KE")
    y = series(rows, "thermo_max_atom_ke_eV")
    if np.any(np.isfinite(y)):
        axs[0, 0].plot(t, y, ls="--", label="thermo max atom KE")
    axs[0, 0].set_ylabel("max atom KE / eV")
    axs[0, 0].legend(fontsize=8)

    axs[0, 1].plot(t, series(rows, "max_vz_A_ps"), marker="o", ms=3, label="dump max vz")
    y = series(rows, "thermo_max_vz_A_ps")
    if np.any(np.isfinite(y)):
        axs[0, 1].plot(t, y, ls="--", label="thermo max vz")
    axs[0, 1].set_ylabel("max vz / A ps^-1")
    axs[0, 1].legend(fontsize=8)

    axs[1, 0].plot(t, series(rows, "surface_zmax_nm"), marker="o", ms=3, label="surface zmax from dump")
    axs[1, 0].plot(t, series(rows, "top1pct_ke_zmin_nm"), marker=".", ms=3, label="top 1% KE z min")
    axs[1, 0].plot(t, series(rows, "top1pct_ke_zmax_nm"), marker=".", ms=3, label="top 1% KE z max")
    axs[1, 0].set_ylabel("z from initial surface / nm")
    axs[1, 0].legend(fontsize=8)

    axs[1, 1].plot(t, series(rows, "thermo_Tl_K"), marker="o", ms=3, label="global thermo T")
    axs[1, 1].plot(t, series(rows, "atom_T_thermal_mean_K"), marker="o", ms=3, label="heated mean thermal T")
    axs[1, 1].plot(t, series(rows, "atom_T_thermal_p95_K"), marker="o", ms=3, label="heated p95 thermal T")
    axs[1, 1].plot(t, series(rows, "atom_T_thermal_max_K"), marker=".", ms=3, label="heated max hotspot ref")
    axs[1, 1].axhline(SI_MELT_K, color="k", ls="--", lw=0.8)
    axs[1, 1].set_ylabel("temperature / K")
    axs[1, 1].legend(fontsize=8)

    for ax in axs.ravel():
        ax.grid(True, alpha=0.25)
        ax.set_xlabel("time / ps")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def plot_depth_temperature_profiles(path: Path, case: str, fields: dict[int, dict[str, object]], targets: list[int]) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 5.2), constrained_layout=True)
    fig.suptitle(f"{case_scope_label(case)}\ndepth-resolved thermal lattice/atom temperature: {case}", fontsize=13)
    plotted = False
    for target in targets:
        item = fields.get(target)
        if item is None:
            continue
        temp = np.asarray(item["temperature"], dtype=float)
        z_edges = np.asarray(item["z_edges"], dtype=float)
        z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
        profile = np.nanmean(temp, axis=1)
        if np.any(np.isfinite(profile)):
            ax.plot(profile, z_centers, label=f"{item['actual_time_ps']:.1f} ps")
            plotted = True
    if not plotted:
        plt.close(fig)
        return
    ax.axvline(SI_MELT_K, color="k", ls="--", lw=0.9, label="Si melt")
    ax.axvline(SI_VAP_K, color="0.35", ls=":", lw=0.9, label="vapor reference")
    ax.axhline(0.0, color="0.4", ls="--", lw=0.8)
    ax.set_xlabel("thermal atom temperature / K")
    ax.set_ylabel("z from initial surface / nm")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def plot_overview(path: Path, case: str, rows: list[dict[str, str]]) -> None:
    t = series(rows, "time_ps")
    fig, axs = plt.subplots(3, 3, figsize=(17, 12), constrained_layout=True)
    fig.suptitle(f"{case_scope_label(case)}\nV8 physical diagnostics overview: {case}", fontsize=16)
    axs[0, 0].plot(t, series(rows, "Te_p95_K"), label="Te p95")
    axs[0, 0].plot(t, series(rows, "thermo_Tl_K"), label="LAMMPS Tl")
    axs[0, 0].plot(t, series(rows, "atom_T_raw_p95_K"), label="atom T raw p95")
    axs[0, 0].plot(t, series(rows, "atom_T_thermal_p95_K"), label="atom T thermal p95")
    axs[0, 0].axhline(SI_MELT_K, color="k", ls="--", lw=0.8)
    axs[0, 0].axhline(SI_VAP_K, color="0.4", ls=":", lw=0.8)
    axs[0, 0].set_ylabel("K")
    axs[0, 0].legend(fontsize=8)
    axs[0, 1].plot(t, series(rows, "low_density_area_fraction"), label="rho<0.7")
    axs[0, 1].plot(t, series(rows, "void_area_fraction"), label="rho<0.3")
    axs[0, 1].set_ylabel("area fraction")
    axs[0, 1].legend(fontsize=8)
    axs[0, 2].plot(t, series(rows, "raw_lift_outlier_sensitive_nm"), label="raw max-z")
    axs[0, 2].plot(t, series(rows, "robust_p95_lift_nm"), label="robust p95")
    axs[0, 2].plot(t, series(rows, "center_crater_depth_nm"), label="center crater")
    axs[0, 2].set_ylabel("nm")
    axs[0, 2].legend(fontsize=8)
    axs[1, 0].plot(t, series(rows, "vapor_candidate_count"), label="vapor candidate - algorithmic only")
    axs[1, 0].plot(t, series(rows, "detached_candidate_count"), label="detached candidate - algorithmic only")
    axs[1, 0].plot(t, series(rows, "stable_ejecta_count"), label="stable ejecta - algorithmic only")
    axs[1, 0].set_ylabel("atom count")
    axs[1, 0].legend(fontsize=8)
    axs[1, 1].plot(t, series(rows, "upward_atoms"), label="upward surface atoms")
    axs[1, 1].plot(t, series(rows, "detached_atom_count"), label="detached by connectivity")
    axs[1, 1].legend(fontsize=8)
    axs[1, 2].plot(t, series(rows, "bulk_component_fraction"), label="largest connected component fraction")
    axs[1, 2].plot(t, series(rows, "ejected_mass_fraction"), label="stable-ejecta mass fraction")
    axs[1, 2].set_ylim(-0.03, 1.03)
    axs[1, 2].set_ylabel("fraction")
    axs[1, 2].legend(fontsize=8)
    axs[2, 0].plot(t, series(rows, "pressure_GPa_p99"), label="P p99")
    axs[2, 0].plot(t, series(rows, "pressure_GPa_p01"), label="P p01")
    axs[2, 0].set_ylabel("GPa")
    axs[2, 0].legend(fontsize=8)
    axs[2, 1].plot(t, series(rows, "von_mises_stress_GPa_p99"), label="von Mises p99")
    axs[2, 1].set_ylabel("GPa")
    axs[2, 1].legend(fontsize=8)
    e = series(rows, "thermo_total_eV")
    finite = np.where(np.isfinite(e))[0]
    if finite.size:
        axs[2, 2].plot(t, e - e[finite[0]], label="Delta total")
    axs[2, 2].plot(t, series(rows, "thermo_ke_eV"), label="atom KE")
    axs[2, 2].plot(t, series(rows, "thermo_pe_eV"), label="atom PE")
    axs[2, 2].legend(fontsize=8)
    for ax in axs.ravel():
        ax.grid(True, alpha=0.25)
        ax.set_xlabel("time / ps")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def plot_surface_profiles(path: Path, case: str, rows: list[dict[str, str]], targets: list[int]) -> None:
    by_time: dict[float, list[tuple[float, float]]] = {}
    for row in rows:
        try:
            time = float(row["time_ps"])
            x = float(row["x_nm"])
            z = float(row["surface_z_nm"])
        except ValueError:
            continue
        nearest = min(targets, key=lambda tt: abs(time - tt))
        if abs(time - nearest) <= 0.6:
            by_time.setdefault(float(nearest), []).append((x, z))
    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)
    for target in sorted(by_time):
        pts = sorted(by_time[target])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], label=f"{target:.0f} ps")
    ax.axhline(0.0, color="k", ls="--", lw=0.8)
    ax.set_title(f"{case_scope_label(case)}\nV8 robust surface profile evolution: {case}")
    ax.set_xlabel("center x / nm")
    ax.set_ylabel("robust surface z / nm")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"[OK] {path}")


def write_readme(path: Path, case: str, rows: list[dict[str, str]], case_dir: Path) -> None:
    max_raw = float(np.nanmax(series(rows, "raw_lift_outlier_sensitive_nm")))
    max_robust = float(np.nanmax(series(rows, "robust_p95_lift_nm")))
    max_stable = int(np.nanmax(series(rows, "stable_ejecta_count")))
    max_detached = int(np.nanmax(series(rows, "detached_candidate_count")))
    max_void = float(np.nanmax(series(rows, "void_area_fraction")))
    max_low = float(np.nanmax(series(rows, "low_density_area_fraction")))
    max_p99 = float(np.nanmax(series(rows, "pressure_GPa_p99")))
    max_surface_t = float(np.nanmax(series(rows, "atom_T_thermal_p95_K")))
    max_global_t = float(np.nanmax(series(rows, "thermo_Tl_K")))
    max_molten_frac = float(np.nanmax(series(rows, "molten_atom_fraction")))
    max_vapor_frac = float(np.nanmax(series(rows, "vapor_threshold_atom_fraction")))
    max_ke = float(np.nanmax(series(rows, "max_atom_ke_eV")))
    max_vz = float(np.nanmax(series(rows, "max_vz_A_ps")))
    max_z = float(np.nanmax(series(rows, "surface_zmax_nm")))
    min_surface_atoms = int(np.nanmin(series(rows, "heated_surface_atom_count")))
    max_surface_atoms = int(np.nanmax(series(rows, "heated_surface_atom_count")))
    init_density = series(rows, "bin_density_mean_check_initial")
    init_density_val = float(init_density[np.isfinite(init_density)][0]) if np.any(np.isfinite(init_density)) else math.nan
    warning = "YES" if max_raw > max_robust * 4.0 + 0.5 else "limited"
    if "ttm_periodic_slab_diagnostic" in case:
        case_scope = (
            "This is a TTM-MD periodic slab diagnostic, not strict free-surface "
            "ablation.  Use it only to verify electron-lattice coupling and "
            "temperature evolution."
        )
    elif "free_surface_ablation" in case:
        case_scope = (
            "This is the free-surface ablation morphology case.  It uses a "
            "non-periodic z boundary and an equivalent near-surface lattice heat "
            "source instead of fix ttm/mod."
        )
    else:
        case_scope = "This is a V8 physical diagnostics case."
    text = f"""# README V8 physical model diagnostics: {case}

## Case scope

{case_scope}

This V8 diagnostic pass is generated from the existing trajectory.  It fixes
the V7 postprocessing definitions before changing the simulation itself.

## V7 issues addressed

- Density now uses `bin_density = current x-z bin count / initial x-z bin count`
  for continuous fields.  Neighbor density is kept only for connectivity and
  ejecta classification.
- Atom temperature now reports both raw kinetic temperature and
  local-drift-subtracted thermal temperature.  Main figures use thermal
  temperature.
- Pressure and von Mises stress use bin-summed LAMMPS `stress/atom` divided by
  bin volume.  LAMMPS documents `compute stress/atom` as pressure*volume, so
  division by a meaningful volume is required before interpreting it as stress.
- Surface lift uses robust p95 profiles after excluding detached and
  low-neighbor-density atoms.

## Current trajectory diagnosis

- Heated surface region definition: atoms from `{HEATED_SURFACE_DEPTH_NM:.3g} nm`
  below the initial surface to `{HEATED_SURFACE_ABOVE_NM:.3g} nm` above it,
  after excluding detached/low-density candidates.  Atom count range in this
  region: `{min_surface_atoms}` to `{max_surface_atoms}`.
- Continuous temperature/density bins require at least
  `{MIN_FIELD_BIN_COUNT}` atoms; stress/pressure bins require at least
  `{MIN_STRESS_BIN_COUNT}` atoms.
- Initial mean bin density check: `{init_density_val:.3g}`; this should be close
  to 1 for the density normalization to be usable.
- Maximum low-density area fraction (`rho/rho0 < 0.7`): `{max_low:.3g}`.
- Maximum void area fraction (`rho/rho0 < 0.3`): `{max_void:.3g}`.
- Maximum raw max-z lift: `{max_raw:.3g} nm`.
- Maximum robust p95 lift: `{max_robust:.3g} nm`.
- Raw-lift outlier warning: `{warning}`.
- Maximum detached candidates: `{max_detached}`.
- Maximum stable ejecta: `{max_stable}`.
- Maximum pressure p99 after bin conversion: `{max_p99:.3g} GPa`.
- Maximum global lattice/atom temperature from thermo: `{max_global_t:.3g} K`.
- Maximum heated-surface thermal temperature p95: `{max_surface_t:.3g} K`.
- Maximum heated-surface thermal temperature max is retained only as a hotspot
  reference.  Do not use the max curve as the primary physical conclusion.
- Maximum molten atom fraction by thermal temperature: `{max_molten_frac:.3g}`.
- Maximum vapor-reference atom fraction by thermal temperature: `{max_vapor_frac:.3g}`.
- Maximum single-atom KE in saved dumps: `{max_ke:.3g} eV`.
- Maximum upward velocity in saved dumps: `{max_vz:.3g} A/ps`.
- Maximum z position relative to the initial surface: `{max_z:.3g} nm`.

## Lost-atoms / heat-source stability checks

This postprocessor also writes:

- `heat_group_diagnostics_v8.csv`: per heat layer atom count, mean/p95
  drift-corrected thermal temperature, max atom KE, max upward velocity, and
  z-range.
- `stability_diagnostics_v8.png`: max KE, max upward velocity, top-1-percent
  KE z-range, surface zmax, pressure/temperature trends before any lost-atoms
  event.

If max KE or max upward velocity jumps sharply around 10-20 ps, treat that as
a heat-source/timestep stability problem first; do not explain it away as only
a top-boundary issue.

## TTM heat/source distribution check

{ttm_source_notes(case_dir, case)}

## Interpretation

Interpret this case according to its declared scope above.  If stable ejecta is
small while raw max-z lift is much larger than robust p95 lift, the apparent
surface motion is dominated by a few high-z atoms rather than a coherent crater
or plume.

For the periodic-slab diagnostic case, do not use these figures to claim a
real crater, stable ejecta, plume, or final morphology.  Its only valid role is
checking whether the TTM-MD coupling, electron temperature output, early lattice
heating, and melt-threshold crossing behave plausibly.

Vapor, detached and stable-ejecta labels in the overview are algorithmic
candidates.  For a periodic slab diagnostic, they are not physical plume or
ejecta conclusions.

## Pressure/stress status

Pressure/stress are now converted using LAMMPS metal-unit virial stress divided
by x-z-bin volume and expressed in GPa.  They are more physical than V7's single
atom-volume proxy, but still require caution because the bin volume is a
coarse-grained control volume and low-count bins are masked.

## What V8 simulation must change next

- The laser source must distinguish incident and absorbed fluence using Si
  reflectivity and optical penetration depth.
- A 30 ps FWHM pulse must be used; 100 ps is an observation window, not a pulse
  duration.
- Lateral size must increase beyond the debug ~4 nm width.
- Run time must increase to 300 ps before judging stable ejecta and plume
  formation.
- Fluence scan should start with 0.5, 0.8 and 1.0 J/cm2, then expand only if the
  simulations remain stable.
"""
    path.write_text(text, encoding="utf-8")
    print(f"[OK] {path}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=Path.cwd())
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="fluence_0p5_100ps")
    ap.add_argument("--equilibration-ps", type=float, default=-1.0, help="Use -1 to infer 2 ps for V61 legacy runs and 0 ps for V80 reset-timestep runs.")
    ap.add_argument("--mass-amu", type=float, default=28.0855)
    args = ap.parse_args()
    project_root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else project_root
    analyze(project_root, output_root, args.case, args.equilibration_ps, args.mass_amu)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
