"""
export_v3b_driver_from_v26.py
=============================

Export a machine-readable V3B driver package from V2.6 threshold metrics
for the main case ep_5p0uj.  Read-only; does not modify metrics CSVs.

Outputs (under results/v26_30ps_threshold_ablation/v3b_driver/):
  v3b_driver_ep5p0uj.json
  v3b_driver_ep5p0uj.csv
  v3b_driver_ep5p0uj_readme.md

Run from project root (after extract_v26):
    .\\.venv\\Scripts\\python.exe scripts\\export_v3b_driver_from_v26.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from extract_v26_30ps_threshold_ablation import MAIN_CASE  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

DRIVER_CASE = MAIN_CASE


def read_metrics_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise RuntimeError(f"empty metrics CSV: {path}")
    out: dict[str, list] = {k: [] for k in rows[0]}
    for row in rows:
        for k, v in row.items():
            if k in ("melt_exists", "vapor_exists"):
                out[k].append(v.strip().lower() == "yes")
            else:
                out[k].append(float(v))
    return {k: np.asarray(v) for k, v in out.items()}


def build_driver_package(metrics: dict[str, np.ndarray]) -> dict:
    k_tmax = int(np.argmax(metrics["Tmax_K"]))
    k_maxd = int(np.argmax(metrics["crater_depth_um"]))

    t_tmax_ps = float(metrics["time_ps"][k_tmax])
    t_maxd_ps = float(metrics["time_ps"][k_maxd])

    driver: dict = {
        "case_name": DRIVER_CASE,
        "source_version": "V2.6",
        "source_driver": "V1.7 local mesh",
        "note_main_case": DRIVER_CASE,
        "t_tmax_peak_ps": round(t_tmax_ps, 3),
        "t_max_crater_depth_ps": round(t_maxd_ps, 3),
        # Tmax peak frame
        "Tmax_peak_K": round(float(metrics["Tmax_K"][k_tmax]), 2),
        "crater_depth_at_tmax_peak_um": round(float(metrics["crater_depth_um"][k_tmax]), 4),
        "crater_radius_at_tmax_peak_um": round(float(metrics["crater_radius_um"][k_tmax]), 4),
        "melt_depth_at_tmax_peak_um": round(float(metrics["melt_depth_um"][k_tmax]), 4),
        "melt_radius_at_tmax_peak_um": round(float(metrics["melt_radius_um"][k_tmax]), 4),
        "vapor_depth_at_tmax_peak_um": round(float(metrics["vapor_depth_um"][k_tmax]), 4),
        "vapor_radius_at_tmax_peak_um": round(float(metrics["vapor_radius_um"][k_tmax]), 4),
        # max crater-depth frame
        "Tmax_at_max_crater_depth_frame_K": round(float(metrics["Tmax_K"][k_maxd]), 2),
        "crater_depth_max_um": round(float(metrics["crater_depth_um"][k_maxd]), 4),
        "crater_radius_at_max_crater_depth_um": round(float(metrics["crater_radius_um"][k_maxd]), 4),
        "melt_depth_at_max_crater_depth_um": round(float(metrics["melt_depth_um"][k_maxd]), 4),
        "melt_radius_at_max_crater_depth_um": round(float(metrics["melt_radius_um"][k_maxd]), 4),
        "vapor_depth_at_max_crater_depth_um": round(float(metrics["vapor_depth_um"][k_maxd]), 4),
        "vapor_radius_at_max_crater_depth_um": round(float(metrics["vapor_radius_um"][k_maxd]), 4),
        # global maxima over time
        "crater_depth_global_max_um": round(float(metrics["crater_depth_um"].max()), 4),
        "crater_radius_global_max_um": round(float(metrics["crater_radius_um"].max()), 4),
        "melt_depth_global_max_um": round(float(metrics["melt_depth_um"].max()), 4),
        "melt_radius_global_max_um": round(float(metrics["melt_radius_um"].max()), 4),
        "vapor_depth_global_max_um": round(float(metrics["vapor_depth_um"].max()), 4),
        "vapor_radius_global_max_um": round(float(metrics["vapor_radius_um"].max()), 4),
        # V3B recommendations
        "recommended_reference_frame_for_geometry": "max_crater_depth_frame",
        "recommended_reference_frame_for_energy": "tmax_peak_frame",
        "recommended_t0_ps": round(t_maxd_ps, 3),
        "recommended_crater_depth_um": round(float(metrics["crater_depth_um"][k_maxd]), 4),
        "recommended_crater_radius_um": round(float(metrics["crater_radius_um"][k_maxd]), 4),
    }
    return driver


def write_driver_readme(path: Path, driver: dict) -> None:
    text = f"""# V3B driver package — {driver['case_name']}

Exported from **V2.6** threshold-equivalent post-processing (source: **V1.7 local mesh**).

## Key reference frames (ep_5p0uj)

| frame | time (ps) | Tmax (K) | crater d (um) | crater r (um) |
| --- | ---: | ---: | ---: | ---: |
| Tmax peak | {driver['t_tmax_peak_ps']} | {driver['Tmax_peak_K']} | {driver['crater_depth_at_tmax_peak_um']} | {driver['crater_radius_at_tmax_peak_um']} |
| max crater depth | {driver['t_max_crater_depth_ps']} | {driver['Tmax_at_max_crater_depth_frame_K']} | {driver['crater_depth_max_um']} | {driver['crater_radius_at_max_crater_depth_um']} |

## V3B usage hints

- **Geometry / crater scale:** use `max_crater_depth_frame` (`recommended_t0_ps` = {driver['recommended_t0_ps']} ps).
- **Thermal / energy peak:** use `tmax_peak_frame` ({driver['t_tmax_peak_ps']} ps, Tmax = {driver['Tmax_peak_K']} K).
- This package is **threshold-equivalent** — not explicit multiphase or ejecta LS-DYNA.

## Files

- `v3b_driver_ep5p0uj.json` — flat key/value driver record
- `v3b_driver_ep5p0uj.csv` — same record, one row
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--v26-root", type=Path,
        default=project_root / "results" / "v26_30ps_threshold_ablation",
    )
    ap.add_argument(
        "--out-dir", type=Path, default=None,
        help="Default: <v26-root>/v3b_driver",
    )
    args = ap.parse_args()

    out_dir = args.out_dir or (args.v26_root / "v3b_driver")
    metrics_path = args.v26_root / DRIVER_CASE / "v26_threshold_metrics.csv"
    if not metrics_path.is_file():
        print(f"[ERROR] metrics not found: {metrics_path}", file=sys.stderr)
        print("        Run scripts\\extract_v26_30ps_threshold_ablation.py first.", file=sys.stderr)
        return 2

    metrics = read_metrics_csv(metrics_path)
    driver = build_driver_package(metrics)
    out_dir.mkdir(parents=True, exist_ok=True)

    json_path = out_dir / "v3b_driver_ep5p0uj.json"
    csv_path = out_dir / "v3b_driver_ep5p0uj.csv"
    readme_path = out_dir / "v3b_driver_ep5p0uj_readme.md"

    json_path.write_text(json.dumps(driver, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(driver.keys()))
        wr.writeheader()
        wr.writerow(driver)
    write_driver_readme(readme_path, driver)

    print("=" * 72)
    print("V3B driver export (from V2.6)")
    print(f"  case     : {DRIVER_CASE}")
    print(f"  out dir  : {out_dir}")
    print(f"  Tmax peak: t={driver['t_tmax_peak_ps']} ps, T={driver['Tmax_peak_K']} K, "
          f"d={driver['crater_depth_at_tmax_peak_um']} um")
    print(f"  max d    : t={driver['t_max_crater_depth_ps']} ps, "
          f"d={driver['crater_depth_max_um']} um, r={driver['crater_radius_at_max_crater_depth_um']} um")
    print("=" * 72)
    print(f"[OK] {json_path}")
    print(f"[OK] {csv_path}")
    print(f"[OK] {readme_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
