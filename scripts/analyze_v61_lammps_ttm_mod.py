"""
analyze_v61_lammps_ttm_mod.py
=============================

Analyze V6B LAMMPS ttm/mod ablation scout trajectories.

Run from project root:
    python scripts\\analyze_v61_lammps_ttm_mod.py --case fluence_0p5_100ps
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

import numpy as np


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


KB_EV_K = 8.617333262145e-5
SI_MELT_K = 1687.0
SI_VAP_K = 3538.0


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
    out: dict[int, dict[str, float]] = {}
    if not log_path.exists():
        return out
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
                out[int(row["Step"])] = row
    return out


def read_te_file(path: Path) -> tuple[float, float, float]:
    vals = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                vals.append(float(parts[3]))
    arr = np.asarray(vals, dtype=float)
    return float(np.nanmin(arr)), float(np.nanmean(arr)), float(np.nanmax(arr))


def te_by_step(case_dir: Path) -> dict[int, tuple[float, float, float]]:
    out = {}
    for p in case_dir.glob("Te_out.*"):
        m = re.search(r"Te_out\.(\d+)$", p.name)
        if m:
            out[int(m.group(1))] = read_te_file(p)
    return out


def nearest_te(step_to_te: dict[int, tuple[float, float, float]], step: int) -> tuple[float, float, float]:
    if not step_to_te:
        return math.nan, math.nan, math.nan
    if step in step_to_te:
        return step_to_te[step]
    key = min(step_to_te.keys(), key=lambda k: abs(k - step))
    return step_to_te[key]


def analyze_case(project_root: Path, output_root: Path, case: str) -> None:
    case_dir = project_root / "results" / "v61_lammps_ttm_mod_ablation_scan" / case
    dump_path = case_dir / f"dump.v61_{case}.lammpstrj"
    log_path = case_dir / "log.lammps"
    if not dump_path.exists():
        raise FileNotFoundError(dump_path)

    thermo = parse_thermo(log_path)
    te = te_by_step(case_dir)
    rows: list[dict[str, str]] = []
    depth_rows: list[dict[str, str]] = []
    initial_surface_A = None

    for step, bounds, cols, data in iter_lammpstrj(dump_path):
        idx = {name: i for i, name in enumerate(cols)}
        z = data[:, idx["z"]]
        vz = data[:, idx["vz"]]
        ke = data[:, idx["c_keatom"]]
        sxx = data[:, idx["c_satom[1]"]]
        syy = data[:, idx["c_satom[2]"]]
        szz = data[:, idx["c_satom[3]"]]

        zmax = float(np.nanmax(z))
        if initial_surface_A is None:
            initial_surface_A = zmax
        lift_nm = (zmax - initial_surface_A) * 0.1
        ejecta_mask = (z > initial_surface_A + 10.0) & (vz > 0.0)
        mobile_mask = ke > 1.0e-12
        temp_atom = np.full_like(ke, np.nan)
        temp_atom[mobile_mask] = 2.0 * ke[mobile_mask] / (3.0 * KB_EV_K)
        stress_proxy_eV = -((sxx + syy + szz) / 3.0)

        th = thermo.get(step, {})
        time_ps = th.get("Time", math.nan)
        thermo_temp = th.get("Temp", math.nan)
        log_press = th.get("Press", math.nan)
        te_min, te_mean, te_max = nearest_te(te, step)

        rows.append({
            "case": case,
            "step": str(step),
            "time_ps": f"{time_ps:.8g}" if math.isfinite(time_ps) else "",
            "natoms": str(len(data)),
            "thermo_temp_K": f"{thermo_temp:.8g}" if math.isfinite(thermo_temp) else "",
            "mean_atom_temperature_K": f"{np.nanmean(temp_atom):.8g}",
            "max_atom_temperature_K": f"{np.nanmax(temp_atom):.8g}",
            "melt_like_atom_count": str(int(np.count_nonzero(temp_atom > SI_MELT_K))),
            "vap_like_atom_count": str(int(np.count_nonzero(temp_atom > SI_VAP_K))),
            "surface_z_A": f"{zmax:.8g}",
            "surface_lift_nm": f"{lift_nm:.8g}",
            "ejecta_atoms": str(int(np.count_nonzero(ejecta_mask))),
            "plume_front_height_nm": f"{max(0.0, lift_nm):.8g}",
            "max_stress_proxy_eV": f"{np.nanmax(stress_proxy_eV):.8g}",
            "min_stress_proxy_eV": f"{np.nanmin(stress_proxy_eV):.8g}",
            "log_pressure": f"{log_press:.8g}" if math.isfinite(log_press) else "",
            "max_electron_temperature_K": f"{te_max:.8g}" if math.isfinite(te_max) else "",
            "mean_electron_temperature_K": f"{te_mean:.8g}" if math.isfinite(te_mean) else "",
            "min_electron_temperature_K": f"{te_min:.8g}" if math.isfinite(te_min) else "",
        })

        bins = np.linspace(0.0, initial_surface_A, 51)
        bin_id = np.digitize(z, bins) - 1
        for b in range(len(bins) - 1):
            mask = bin_id == b
            if not np.any(mask):
                continue
            temp_mean = np.nanmean(temp_atom[mask])
            depth_mid_nm = (initial_surface_A - 0.5 * (bins[b] + bins[b + 1])) * 0.1
            depth_rows.append({
                "case": case,
                "step": str(step),
                "time_ps": f"{time_ps:.8g}" if math.isfinite(time_ps) else "",
                "depth_mid_nm": f"{depth_mid_nm:.8g}",
                "atom_count": str(int(np.count_nonzero(mask))),
                "mean_atom_temperature_K": f"{temp_mean:.8g}" if math.isfinite(temp_mean) else "",
                "mean_stress_proxy_eV": f"{np.nanmean(stress_proxy_eV[mask]):.8g}",
            })

    out_dir = output_root / "lightweight_results" / "v61_lammps_ttm_mod_ablation_scan"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / f"v61_{case}_summary.csv"
    depth_path = out_dir / f"v61_{case}_depth_profiles.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with depth_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(depth_rows[0].keys()))
        writer.writeheader()
        writer.writerows(depth_rows)
    print(f"[OK] {summary_path}")
    print(f"[OK] {depth_path}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=Path.cwd())
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="fluence_0p5_100ps")
    args = ap.parse_args()
    project_root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else project_root
    analyze_case(project_root, output_root, args.case)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
