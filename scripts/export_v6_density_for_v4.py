"""
export_v6_density_for_v4.py
===========================

Bin V6 LAMMPS snapshot atoms into a 2D density grid for V4 structured-light handoff.

Output:
  results/v6_lammps_ttm_md_pilot/<case>/v6_density_grid_<t>ps.csv
  results/v6_lammps_ttm_md_pilot/<case>/v6_v4_handoff_manifest.json

Run after LAMMPS snapshots exist:
    python scripts\\export_v6_density_for_v4.py --case ep_5p0uj
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import tomllib
from pathlib import Path

_THIS = Path(__file__).resolve().parent
if str(_THIS) not in sys.path:
    sys.path.insert(0, str(_THIS))

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

from plot_v6_lammps_ttm_snapshots import find_snapshots, read_lammpstrj, _view_coords  # noqa: E402


def bin_density(r_um, z_um, r_edges, z_edges) -> np.ndarray:
    H, _, _ = np.histogram2d(r_um, z_um, bins=[r_edges, z_edges])
    return H.T


def export_case(case: str, root: Path, cfg: dict, *, laser_axis: str = "+x") -> int:
    res_dir = root / "results/v6_lammps_ttm_md_pilot" / case
    case_dir = root / "models/v6_lammps_ttm_md_pilot" / case
    snaps = find_snapshots(res_dir)
    if not snaps:
        snaps = find_snapshots(case_dir)
    if not snaps:
        print(f"[WARN] no snapshots in {res_dir}", file=sys.stderr)
        return 2

    exp = cfg["export"]
    r_min, r_max = float(exp["r_min_um"]), float(exp["r_max_um"])
    z_min, z_max = float(exp["z_min_um"]), float(exp["z_max_um"])
    dr = float(exp["bin_r_um"])
    dz = float(exp["bin_z_um"])
    r_edges = np.arange(r_min, r_max + dr, dr)
    z_edges = np.arange(z_min, z_max + dz, dz)
    r_cent = 0.5 * (r_edges[:-1] + r_edges[1:])
    z_cent = 0.5 * (z_edges[:-1] + z_edges[1:])

    manifest = {"case": case, "frames": []}
    for t_ps, path in snaps:
        x, y, z = read_lammpstrj(path)
        r_um, z_um, _ = _view_coords(x, y, z, laser_axis)
        H = bin_density(r_um, z_um, r_edges, z_edges)
        out_csv = res_dir / f"v6_density_grid_{int(t_ps)}ps.csv"
        with out_csv.open("w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(["r_um", "z_um", "number_density_proxy"])
            for iz, zz in enumerate(z_cent):
                for ir, rr in enumerate(r_cent):
                    if H[iz, ir] > 0:
                        wr.writerow([f"{rr:.4f}", f"{zz:.4f}", f"{H[iz, ir]:.6f}"])
        manifest["frames"].append({
            "time_ps": t_ps,
            "csv": str(out_csv.relative_to(root)).replace("\\", "/"),
            "n_atoms": int(len(x)),
        })
        print(f"[OK] {out_csv}  atoms={len(x)}")

    man_path = res_dir / "v6_v4_handoff_manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] {man_path}")
    return 0


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", default="ep_5p0uj")
    ap.add_argument("--toml", type=Path, default=root / "config/v6_lammps_ttm_md_pilot.toml")
    args = ap.parse_args()

    case_dir = root / "models/v6_lammps_ttm_md_pilot" / args.case
    laser_axis = "+x"
    manifest = case_dir / "v6_run_manifest.json"
    if manifest.is_file():
        laser_axis = str(json.loads(manifest.read_text(encoding="utf-8")).get("laser_axis", laser_axis))

    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)
    return export_case(args.case, root, cfg, laser_axis=laser_axis)


if __name__ == "__main__":
    raise SystemExit(main())
