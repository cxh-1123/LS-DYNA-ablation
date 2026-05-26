"""
check_v17_outputs.py
====================

Inspect V1.7 LS-DYNA results for selected 30 ps cases.

Checks:
  * Normal termination in messag
  * tprint parses
  * Tmax(t), melt/vapor thresholds
  * peak location near r=0, top surface (axis_top_node_id from registry)

Output:
  results/v17_30ps_local/v17_case_summary.csv

Run:
    .\\.venv\\Scripts\\python.exe scripts\\check_v17_outputs.py
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import tomllib
from pathlib import Path

import numpy as np

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from lsdyna_tprint import parse_tprint  # noqa: E402


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


REQUIRED = ["d3plot", "d3hsp", "messag", "tprint"]


def messag_ok(path: Path) -> bool:
    raw = path.read_text(encoding="utf-8", errors="ignore")
    return "normaltermination" in re.sub(r"\s+", "", raw).lower()


def load_mesh_nodes(path: Path) -> dict[int, tuple[int, int]]:
    """node_id -> (i, j) grid indices."""
    mapping: dict[int, tuple[int, int]] = {}
    with path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            mapping[int(row["node_id"])] = (int(row["i"]), int(row["j"]))
    return mapping


def tprint_to_grid(T_vec: np.ndarray, node_map: dict[int, tuple[int, int]],
                   NR: int, NZ: int) -> np.ndarray:
    """Map 1-D tprint node vector (0-indexed by node_id-1) to (NZ+1, NR+1)."""
    grid = np.full((NZ + 1, NR + 1), np.nan)
    for nid, (i, j) in node_map.items():
        idx = nid - 1
        if 0 <= idx < len(T_vec):
            grid[j, i] = T_vec[idx]
    return grid


def check_case(row: dict, project_root: Path, T_melt: float, T_vap: float,
               node_map: dict[int, tuple[int, int]]) -> dict:
    name = row["name"]
    NR = int(row["NR"])
    NZ = int(row["NZ"])
    axis_top = int(row["axis_top_node_id"])
    run_flag = row["run_lsdyna_final"].strip().lower() in ("true", "1", "yes")

    out = {
        "case": name,
        "Ep_uJ": row["Ep_uJ"],
        "tau_ps": row["tau_ps"],
        "ran_lsdyna": "no",
        "status": "",
        "normal_termination": "no",
        "T_peak_K": "",
        "t_peak_ns": "",
        "node_at_peak": "",
        "peak_at_axis_top": "",
        "melt_reached": "",
        "vapor_reached": "",
        "observed_regime": "",
        "notes": row.get("notes", ""),
    }

    if not run_flag:
        out["status"] = "not scheduled"
        return out

    out_dir = project_root / "results" / "v17_30ps_local" / name
    if not out_dir.is_dir():
        out["status"] = f"missing result dir: {out_dir.relative_to(project_root)}"
        return out

    out["ran_lsdyna"] = "yes"
    missing = [f for f in REQUIRED if not (out_dir / f).is_file()]
    if missing:
        out["status"] = f"missing: {', '.join(missing)}"
        return out

    if not messag_ok(out_dir / "messag"):
        out["status"] = "no Normal termination"
        return out
    out["normal_termination"] = "yes"

    try:
        times_ms, T_all = parse_tprint(out_dir / "tprint")
    except Exception as e:
        out["status"] = f"tprint parse failed: {e}"
        return out

    flat_idx = int(np.nanargmax(T_all))
    n_nodes = T_all.shape[1]
    snap_k = flat_idx // n_nodes
    node_idx = flat_idx % n_nodes
    node_id = node_idx + 1
    T_peak = float(T_all[snap_k, node_idx])
    t_peak_ns = float(times_ms[snap_k]) * 1e6

    grid = tprint_to_grid(T_all[snap_k], node_map, NR, NZ)
    i_pk, j_pk = node_map.get(node_id, (-1, -1))
    peak_at_top = (node_id == axis_top) or (i_pk == 0 and j_pk == NZ)

    out["T_peak_K"] = f"{T_peak:.2f}"
    out["t_peak_ns"] = f"{t_peak_ns:.4f}"
    out["node_at_peak"] = str(node_id)
    out["peak_at_axis_top"] = "yes" if peak_at_top else "no"
    out["melt_reached"] = "yes" if T_peak >= T_melt else "no"
    out["vapor_reached"] = "yes" if T_peak >= T_vap else "no"
    if T_peak >= T_vap:
        out["observed_regime"] = "vapor/ablation candidate"
    elif T_peak >= T_melt:
        out["observed_regime"] = "melt only"
    else:
        out["observed_regime"] = "no melt"
    out["status"] = "OK"
    return out


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--toml", type=Path,
                    default=project_root / "config" / "v17_30ps_local_mesh.toml")
    ap.add_argument("--registry", type=Path,
                    default=project_root / "models" / "v17_30ps_local" / "v17_case_registry.csv")
    ap.add_argument("--mesh-nodes", type=Path,
                    default=project_root / "models" / "v17_30ps_local" / "v17_mesh_nodes.csv")
    ap.add_argument("--out-root", type=Path,
                    default=project_root / "results" / "v17_30ps_local")
    args = ap.parse_args()

    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)
    T_melt = float(cfg["material"]["T_melt_K"])
    T_vap = float(cfg["material"]["T_vap_K"])

    if not args.registry.is_file():
        print(f"[ERROR] registry missing: {args.registry}", file=sys.stderr)
        return 2
    if not args.mesh_nodes.is_file():
        print(f"[ERROR] mesh nodes missing: {args.mesh_nodes}", file=sys.stderr)
        return 3

    node_map = load_mesh_nodes(args.mesh_nodes)
    with args.registry.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    summary = [check_case(r, project_root, T_melt, T_vap, node_map) for r in rows]
    args.out_root.mkdir(parents=True, exist_ok=True)
    out_csv = args.out_root / "v17_case_summary.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        wr.writeheader()
        wr.writerows(summary)

    print("=" * 78)
    print("V1.7 output check")
    print(f"  summary -> {out_csv}")
    print("=" * 78)
    for s in summary:
        print(f"  {s['case']:12}  status={s['status']:20}  "
              f"T_peak={s['T_peak_K'] or 'n/a':>8}  "
              f"regime={s['observed_regime'] or 'n/a'}  "
              f"axis_top={s['peak_at_axis_top'] or 'n/a'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
