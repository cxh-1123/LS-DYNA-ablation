"""
analyze_v60_lammps_ttm_md.py
============================

Lightweight analysis for V6A LAMMPS TTM-MD trajectory files.

Run from project root:
    python scripts\\analyze_v60_lammps_ttm_md.py --project-root D:\\cxh-daima\\LS-DYNA-ablation
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


def parse_thermo_times(log_path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
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
                    step = int(float(parts[0]))
                    time = float(parts[1])
                except ValueError:
                    header = None
                    continue
                out[step] = time
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
        if not m:
            continue
        out[int(m.group(1))] = read_te_file(p)
    return out


def analyze_case(project_root: Path, case: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    case_dir = project_root / "results" / "v60_lammps_ttm_md_pilot" / case
    dump_path = case_dir / f"dump.v60_{case}.lammpstrj"
    log_path = case_dir / "log.lammps"
    if not dump_path.exists():
        raise FileNotFoundError(dump_path)

    step_to_time = parse_thermo_times(log_path)
    step_to_te = te_by_step(case_dir)
    rows = []
    depth_rows = []
    initial_surface_A = None

    for step, bounds, cols, data in iter_lammpstrj(dump_path):
        idx = {name: i for i, name in enumerate(cols)}
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

        zmax = float(np.nanmax(z))
        if initial_surface_A is None:
            initial_surface_A = zmax
        ejecta_mask = (z > initial_surface_A + 10.0) & (vz > 0.0)
        mobile_mask = ke > 1.0e-12
        temp_atom = np.full_like(ke, np.nan)
        temp_atom[mobile_mask] = 2.0 * ke[mobile_mask] / (3.0 * KB_EV_K)

        # LAMMPS stress/atom is a virial quantity with pressure*volume units.
        # Keep V6A as a signed stress proxy; V6B will bin atoms into real local
        # volumes before reporting calibrated pressure in GPa.
        stress_proxy_eV = -((sxx + syy + szz) / 3.0)

        te_min, te_mean, te_max = step_to_te.get(step, (math.nan, math.nan, math.nan))
        time_ps = step_to_time.get(step, math.nan)

        rows.append({
            "case": case,
            "step": str(step),
            "time_ps": f"{time_ps:.8g}" if math.isfinite(time_ps) else "",
            "natoms": str(len(data)),
            "surface_z_A": f"{zmax:.8g}",
            "surface_lift_nm": f"{(zmax - initial_surface_A) * 0.1:.8g}",
            "ejecta_atoms": str(int(np.count_nonzero(ejecta_mask))),
            "plume_front_height_nm": f"{max(0.0, (zmax - initial_surface_A) * 0.1):.8g}",
            "max_atom_temperature_K": f"{np.nanmax(temp_atom):.8g}",
            "mean_atom_temperature_K": f"{np.nanmean(temp_atom):.8g}",
            "max_stress_proxy_eV": f"{np.nanmax(stress_proxy_eV):.8g}",
            "min_stress_proxy_eV": f"{np.nanmin(stress_proxy_eV):.8g}",
            "max_electron_temperature_K": f"{te_max:.8g}" if math.isfinite(te_max) else "",
            "mean_electron_temperature_K": f"{te_mean:.8g}" if math.isfinite(te_mean) else "",
            "min_electron_temperature_K": f"{te_min:.8g}" if math.isfinite(te_min) else "",
        })

        bins = np.linspace(0.0, initial_surface_A, 41)
        bin_id = np.digitize(z, bins) - 1
        for b in range(len(bins) - 1):
            mask = bin_id == b
            if not np.any(mask):
                continue
            depth_mid_nm = (initial_surface_A - 0.5 * (bins[b] + bins[b + 1])) * 0.1
            depth_rows.append({
                "case": case,
                "step": str(step),
                "time_ps": f"{time_ps:.8g}" if math.isfinite(time_ps) else "",
                "depth_mid_nm": f"{depth_mid_nm:.8g}",
                "atom_count": str(int(np.count_nonzero(mask))),
                "mean_atom_temperature_K": f"{np.nanmean(temp_atom[mask]):.8g}",
                "mean_stress_proxy_eV": f"{np.nanmean(stress_proxy_eV[mask]):.8g}",
            })

    return rows, depth_rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=Path.cwd())
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="ttm_initial_pulse_5uj_equiv")
    args = ap.parse_args()

    project_root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else project_root
    out_dir = output_root / "lightweight_results" / "v60_lammps_ttm_md_pilot"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows, depth_rows = analyze_case(project_root, args.case)
    summary_path = out_dir / "v60_ttm_md_summary.csv"
    depth_path = out_dir / "v60_depth_profiles.csv"

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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
