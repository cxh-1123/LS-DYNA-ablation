"""
extract_v26_30ps_threshold_ablation.py
======================================

V2.6 -- temperature-threshold equivalent ablation post-processing for the
V1.7 local refined 30 ps LS-DYNA model.

Reads V1.7 tprint files (non-uniform r-z mesh) and computes melt / vapor /
crater metrics at each output time.  Python-only; does NOT modify .k files,
d3plot, or re-run LS-DYNA.

Outputs:
  results/v26_30ps_threshold_ablation/<case>/v26_threshold_metrics.csv
  results/v26_30ps_threshold_ablation/v26_case_summary.csv

Run from project root:
    .\\.venv\\Scripts\\python.exe scripts\\extract_v26_30ps_threshold_ablation.py
"""

from __future__ import annotations

import argparse
import csv
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from lsdyna_tprint import parse_tprint  # noqa: E402
from check_v17_outputs import load_mesh_nodes, tprint_to_grid  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

MAIN_CASE = "ep_5p0uj"
DEFAULT_CASE_ORDER = ["ep_0p5uj", "ep_1p0uj", "ep_2p0uj", "ep_5p0uj"]

METRICS_FIELDS = [
    "time_ns",
    "time_ps",
    "Tmax_K",
    "melt_exists",
    "vapor_exists",
    "melt_depth_um",
    "melt_radius_um",
    "vapor_depth_um",
    "vapor_radius_um",
    "crater_depth_um",
    "crater_radius_um",
    "r_at_Tmax_um",
    "z_at_Tmax_um",
]

SUMMARY_FIELDS = [
    "case",
    "status",
    "Tmax_peak_K",
    "t_peak_ns",
    "t_peak_ps",
    "max_melt_depth_um",
    "max_melt_radius_um",
    "max_vapor_depth_um",
    "max_vapor_radius_um",
    "max_crater_depth_um",
    "max_crater_radius_um",
    "final_regime",
    "peak_at_axis_top",
    "peak_warning",
]


@dataclass(frozen=True)
class V17MeshGrid:
    r_um: np.ndarray
    z_um: np.ndarray
    node_map: dict[int, tuple[int, int]]
    NR: int
    NZ: int
    axis_top_node_id: int

    @property
    def z_top_um(self) -> float:
        return float(self.z_um[-1])

    @property
    def n_nodes_expected(self) -> int:
        return (self.NR + 1) * (self.NZ + 1)


def load_mesh_grid(mesh_csv: Path, registry_row: dict) -> V17MeshGrid:
    r_um, z_um = [], []
    node_map = load_mesh_nodes(mesh_csv)
    with mesh_csv.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            i, j = int(row["i"]), int(row["j"])
            if len(r_um) <= i:
                r_um.append(float(row["r_um"]))
            if len(z_um) <= j:
                z_um.append(float(row["z_um"]))
    NR = int(registry_row["NR"])
    NZ = int(registry_row["NZ"])
    return V17MeshGrid(
        r_um=np.asarray(r_um, dtype=float),
        z_um=np.asarray(z_um, dtype=float),
        node_map=node_map,
        NR=NR,
        NZ=NZ,
        axis_top_node_id=int(registry_row["axis_top_node_id"]),
    )


def interp_depth_axis_um(T_col: np.ndarray, z_um: np.ndarray, T_thresh: float) -> float:
    """
    Depth from top surface on r=0 axis where T crosses T_thresh (linear interp).
    """
    j_top = len(z_um) - 1
    z_top = float(z_um[j_top])
    if T_col[j_top] < T_thresh:
        return 0.0
    for j in range(j_top - 1, -1, -1):
        if T_col[j] < T_thresh:
            t_above = T_col[j + 1]
            t_below = T_col[j]
            denom = t_above - t_below
            frac = (t_above - T_thresh) / denom if denom > 0.0 else 0.0
            z_cross = z_um[j + 1] - frac * (z_um[j + 1] - z_um[j])
            return float(z_top - z_cross)
    return float(z_top - z_um[0])


def interp_radius_top_um(T_row: np.ndarray, r_um: np.ndarray, T_thresh: float) -> float:
    """Radius on top surface where T drops below T_thresh (linear interp)."""
    if T_row[0] < T_thresh:
        return 0.0
    for i in range(1, len(r_um)):
        if T_row[i] < T_thresh:
            t_in = T_row[i - 1]
            t_out = T_row[i]
            denom = t_in - t_out
            frac = (t_in - T_thresh) / denom if denom > 0.0 else 0.0
            return float(r_um[i - 1] + frac * (r_um[i] - r_um[i - 1]))
    return float(r_um[-1])


def threshold_zone_metrics(
    T_grid: np.ndarray,
    mesh: V17MeshGrid,
    T_thresh: float,
) -> dict:
    exists = bool((T_grid >= T_thresh).any())
    if not exists:
        return {"depth_um": 0.0, "radius_um": 0.0, "area_um2": 0.0, "exists": False}

    j_top = mesh.NZ
    depth_um = interp_depth_axis_um(T_grid[:, 0], mesh.z_um, T_thresh)
    radius_um = interp_radius_top_um(T_grid[j_top, :], mesh.r_um, T_thresh)

    dr = np.diff(mesh.r_um)
    dz = np.diff(mesh.z_um)
    Tc = 0.25 * (
        T_grid[:-1, :-1] + T_grid[1:, :-1]
        + T_grid[:-1, 1:] + T_grid[1:, 1:]
    )
    mask = Tc >= T_thresh
    area_um2 = float(np.sum(dr[np.newaxis, :] * dz[:, np.newaxis] * mask))
    return {
        "depth_um": float(depth_um),
        "radius_um": float(radius_um),
        "area_um2": area_um2,
        "exists": True,
    }


def peak_location_um(T_grid: np.ndarray, mesh: V17MeshGrid) -> tuple[float, float, int, int]:
    j_pk, i_pk = np.unravel_index(int(np.nanargmax(T_grid)), T_grid.shape)
    return float(mesh.r_um[i_pk]), float(mesh.z_um[j_pk]), int(i_pk), int(j_pk)


def final_regime_from_Tmax(Tmax_peak: float, T_melt: float, T_vap: float) -> str:
    if Tmax_peak >= T_vap:
        return "vapor/ablation candidate"
    if Tmax_peak >= T_melt:
        return "melt only"
    return "no melt"


def process_case(
    case_name: str,
    case_dir: Path,
    mesh: V17MeshGrid,
    out_root: Path,
    T_melt: float,
    T_vap: float,
) -> dict | None:
    tprint = case_dir / "tprint"
    if not case_dir.is_dir():
        print(f"  [SKIP] {case_name}: result dir not found", file=sys.stderr)
        return None
    if not tprint.is_file():
        print(f"  [SKIP] {case_name}: no tprint", file=sys.stderr)
        return None

    try:
        times_ms, T_all = parse_tprint(tprint)
    except Exception as exc:
        print(f"  [FAIL] {case_name}: tprint parse failed: {exc}", file=sys.stderr)
        return None

    if T_all.shape[1] != mesh.n_nodes_expected:
        print(
            f"  [FAIL] {case_name}: node count {T_all.shape[1]} != "
            f"expected {mesh.n_nodes_expected} (mesh_nodes mismatch)",
            file=sys.stderr,
        )
        return None

    times_ns = np.asarray(times_ms, dtype=float) * 1e6
    rows: list[dict] = []
    Tmax_track = -np.inf
    t_peak_ns = 0.0
    peak_axis_ok = True
    peak_warning = ""

    for k in range(len(times_ns)):
        T_grid = tprint_to_grid(T_all[k], mesh.node_map, mesh.NR, mesh.NZ)
        Tmax_K = float(np.nanmax(T_grid))
        if Tmax_K > Tmax_track:
            Tmax_track = Tmax_K
            t_peak_ns = float(times_ns[k])
            r_pk, z_pk, i_pk, j_pk = peak_location_um(T_grid, mesh)
            node_id = i_pk + j_pk * (mesh.NR + 1) + 1
            peak_axis_ok = (
                node_id == mesh.axis_top_node_id
                or (i_pk == 0 and j_pk == mesh.NZ)
            )

        melt = threshold_zone_metrics(T_grid, mesh, T_melt)
        vapor = threshold_zone_metrics(T_grid, mesh, T_vap)
        r_at, z_at, _, _ = peak_location_um(T_grid, mesh)

        rows.append({
            "time_ns": float(times_ns[k]),
            "time_ps": float(times_ns[k] * 1e3),
            "Tmax_K": Tmax_K,
            "melt_exists": "yes" if melt["exists"] else "no",
            "vapor_exists": "yes" if vapor["exists"] else "no",
            "melt_depth_um": melt["depth_um"],
            "melt_radius_um": melt["radius_um"],
            "vapor_depth_um": vapor["depth_um"],
            "vapor_radius_um": vapor["radius_um"],
            "crater_depth_um": vapor["depth_um"],
            "crater_radius_um": vapor["radius_um"],
            "r_at_Tmax_um": r_at,
            "z_at_Tmax_um": z_at,
        })

    if not peak_axis_ok:
        peak_warning = "peak not at r=0 top surface (axis_top)"
        print(f"  [WARN] {case_name}: {peak_warning}", file=sys.stderr)

    case_out = out_root / case_name
    case_out.mkdir(parents=True, exist_ok=True)
    metrics_csv = case_out / "v26_threshold_metrics.csv"
    with metrics_csv.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=METRICS_FIELDS)
        wr.writeheader()
        wr.writerows(rows)
    print(f"  [ OK ] {case_name}: {len(rows)} snapshots -> {metrics_csv}")

    arr_melt_d = np.array([r["melt_depth_um"] for r in rows])
    arr_melt_r = np.array([r["melt_radius_um"] for r in rows])
    arr_vapor_d = np.array([r["vapor_depth_um"] for r in rows])
    arr_vapor_r = np.array([r["vapor_radius_um"] for r in rows])

    regime = final_regime_from_Tmax(Tmax_track, T_melt, T_vap)
    return {
        "case": case_name,
        "status": "OK",
        "Tmax_peak_K": f"{Tmax_track:.2f}",
        "t_peak_ns": f"{t_peak_ns:.6f}",
        "t_peak_ps": f"{t_peak_ns * 1e3:.3f}",
        "max_melt_depth_um": f"{arr_melt_d.max():.4f}",
        "max_melt_radius_um": f"{arr_melt_r.max():.4f}",
        "max_vapor_depth_um": f"{arr_vapor_d.max():.4f}",
        "max_vapor_radius_um": f"{arr_vapor_r.max():.4f}",
        "max_crater_depth_um": f"{arr_vapor_d.max():.4f}",
        "max_crater_radius_um": f"{arr_vapor_r.max():.4f}",
        "final_regime": regime,
        "peak_at_axis_top": "yes" if peak_axis_ok else "no",
        "peak_warning": peak_warning,
    }


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--toml", type=Path,
                    default=project_root / "config" / "v17_30ps_local_mesh.toml")
    ap.add_argument("--registry", type=Path,
                    default=project_root / "models" / "v17_30ps_local" / "v17_case_registry.csv")
    ap.add_argument("--mesh-nodes", type=Path,
                    default=project_root / "models" / "v17_30ps_local" / "v17_mesh_nodes.csv")
    ap.add_argument("--v17-results-root", type=Path,
                    default=project_root / "results" / "v17_30ps_local")
    ap.add_argument("--out-root", type=Path,
                    default=project_root / "results" / "v26_30ps_threshold_ablation")
    ap.add_argument("--cases", nargs="+", default=DEFAULT_CASE_ORDER)
    args = ap.parse_args()

    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)
    T_melt = float(cfg["material"]["T_melt_K"])
    T_vap = float(cfg["material"]["T_vap_K"])
    T_init = float(cfg["material"]["T_init_K"])

    if not args.registry.is_file():
        print(f"[ERROR] registry missing: {args.registry}", file=sys.stderr)
        return 2
    if not args.mesh_nodes.is_file():
        print(f"[ERROR] mesh_nodes missing: {args.mesh_nodes}", file=sys.stderr)
        return 2

    with args.registry.open("r", encoding="utf-8") as fh:
        registry = {row["name"]: row for row in csv.DictReader(fh)}

    print("=" * 78)
    print("V2.6 -- 30 ps threshold-equivalent ablation (driver: V1.7 local mesh)")
    print(f"  main case      : {MAIN_CASE}")
    print(f"  V1.7 results   : {args.v17_results_root}")
    print(f"  output root    : {args.out_root}")
    print(f"  T_init/T_m/T_v : {T_init}/{T_melt}/{T_vap} K")
    print("=" * 78)

    args.out_root.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict] = []
    missing_rows: list[dict] = []

    for case_name in args.cases:
        if case_name not in registry:
            print(f"  [SKIP] {case_name}: not in registry", file=sys.stderr)
            missing_rows.append({
                "case": case_name, "status": "not in registry",
                "final_regime": "missing",
            })
            continue
        mesh = load_mesh_grid(args.mesh_nodes, registry[case_name])
        row = process_case(
            case_name,
            args.v17_results_root / case_name,
            mesh,
            args.out_root,
            T_melt,
            T_vap,
        )
        if row is not None:
            summary_rows.append(row)
        else:
            missing_rows.append({
                "case": case_name,
                "status": "not run",
                "Tmax_peak_K": "",
                "t_peak_ns": "",
                "t_peak_ps": "",
                "max_melt_depth_um": "",
                "max_melt_radius_um": "",
                "max_vapor_depth_um": "",
                "max_vapor_radius_um": "",
                "max_crater_depth_um": "",
                "max_crater_radius_um": "",
                "final_regime": "missing",
                "peak_at_axis_top": "",
                "peak_warning": "",
            })

    if not summary_rows:
        print("\n[ERROR] No cases produced metrics.", file=sys.stderr)
        return 3

    all_summary = summary_rows + missing_rows
    summary_csv = args.out_root / "v26_case_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        wr.writeheader()
        for row in all_summary:
            wr.writerow({k: row.get(k, "") for k in SUMMARY_FIELDS})

    print()
    print("=" * 78)
    print("Per-case maxima")
    print("=" * 78)
    print(f"  {'case':12} {'Tpeak(K)':>10} {'melt_d':>9} {'melt_r':>9} "
          f"{'vap_d':>9} {'vap_r':>9} {'regime'}")
    for row in summary_rows:
        print(
            f"  {row['case']:12} {row['Tmax_peak_K']:>10} "
            f"{row['max_melt_depth_um']:>9} {row['max_melt_radius_um']:>9} "
            f"{row['max_vapor_depth_um']:>9} {row['max_vapor_radius_um']:>9} "
            f"{row['final_regime']}"
        )
    for row in missing_rows:
        print(f"  {row['case']:12}  -- missing / not run")

    print(f"\n[DONE] summary -> {summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
