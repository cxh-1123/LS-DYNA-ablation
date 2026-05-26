"""
build_v17_30ps_local_mesh.py
============================

V1.7 step 1 -- generate local refined 2D axisymmetric LS-DYNA .k files
for the 30 ps thermal model.

Reads:
  - config/v17_30ps_local_mesh.toml

Writes:
  - models/v17_30ps_local/v17_<case>.k
  - models/v17_30ps_local/v17_case_registry.csv
  - models/v17_30ps_local/v17_mesh_summary.csv
  - models/v17_30ps_local/v17_mesh_nodes.csv

Does NOT run LS-DYNA.

Run from project root:
    .\\.venv\\Scripts\\python.exe scripts\\build_v17_30ps_local_mesh.py
"""

from __future__ import annotations

import argparse
import csv
import sys
import tomllib
from pathlib import Path

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from _build_v17_case import (  # noqa: E402
    V17Case, build_k_text, build_mesh, write_k, node_id_at, axis_top_node_id,
    resolve_database_dt_ms,
)


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


def write_mesh_nodes_csv(mesh, path: Path) -> None:
    rows = []
    for j in range(mesh.NZ + 1):
        for i in range(mesh.NR + 1):
            nid = node_id_at(mesh, i, j)
            rows.append({
                "node_id": nid,
                "i": i,
                "j": j,
                "r_um": f"{mesh.r_um[i]:.6f}",
                "z_um": f"{mesh.z_um[j]:.6f}",
                "r_mm": f"{mesh.r_mm[i]:.10e}",
                "z_mm": f"{mesh.z_mm[j]:.10e}",
            })
    with path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--toml", type=Path,
        default=project_root / "config" / "v17_30ps_local_mesh.toml",
    )
    ap.add_argument(
        "--out-dir", type=Path,
        default=project_root / "models" / "v17_30ps_local",
    )
    args = ap.parse_args()

    if not args.toml.is_file():
        print(f"[ERROR] config not found: {args.toml}", file=sys.stderr)
        return 2
    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)

    mesh = build_mesh(cfg)
    st = mesh.stats
    mat = cfg["material"]
    laser = cfg["laser"]
    sol = cfg["solver"]

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("V1.7 -- build local refined 30 ps mesh + .k files")
    print(f"  r = 0 .. {mesh.r_um[-1]:g} um")
    print(f"  z = 0 .. {mesh.z_um[-1]:g} um")
    print(f"  NR x NZ = {st['NR']} x {st['NZ']}")
    print(f"  nodes = {st['n_nodes']},  elements = {st['n_elements']}")
    print(f"  dr = {st['dr_min_um']:.4g} .. {st['dr_max_um']:.4g} um")
    print(f"  dz_min = {st['dz_min_nm']:.4g} nm  (at top surface)")
    print(f"  axis-top node id = {axis_top_node_id(mesh)}")
    t_end_ms = float(sol["t_end_ns"]) * 1e-6
    dt_ms, dt_desc = resolve_database_dt_ms(cfg, t_end_ms)
    print(f"  database DT = {dt_ms * 1e6:g} ns  ({dt_desc})")
    print("=" * 78)

    # mesh summary
    summary_path = args.out_dir / "v17_mesh_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=[
            "radius_um", "thickness_um", "NR", "NZ", "n_nodes", "n_elements",
            "dr_min_um", "dr_max_um", "dz_min_nm", "dz_max_nm",
            "axis_top_node_id",
        ])
        wr.writeheader()
        wr.writerow({
            "radius_um": mesh.r_um[-1],
            "thickness_um": mesh.z_um[-1],
            "NR": st["NR"], "NZ": st["NZ"],
            "n_nodes": st["n_nodes"], "n_elements": st["n_elements"],
            "dr_min_um": f"{st['dr_min_um']:.6f}",
            "dr_max_um": f"{st['dr_max_um']:.6f}",
            "dz_min_nm": f"{st['dz_min_nm']:.4f}",
            "dz_max_nm": f"{st['dz_max_um']*1000:.4f}",
            "axis_top_node_id": axis_top_node_id(mesh),
        })
    write_mesh_nodes_csv(mesh, args.out_dir / "v17_mesh_nodes.csv")
    print(f"[OK] {summary_path}")
    print(f"[OK] {args.out_dir / 'v17_mesh_nodes.csv'}")

    registry_rows = []
    for sc in cfg["selected_cases"]:
        case = V17Case(
            name=sc["name"],
            Ep_uJ=float(sc["Ep_uJ"]),
            run_lsdyna=bool(sc["run_lsdyna_hint"]),
            notes=str(sc.get("notes", "")),
            tau_ps=float(laser["pulse_width_ps"]),
            w0_um=float(laser["spot_radius_um"]),
            A=float(laser["absorption_A"]),
            t_end_ns=float(sol["t_end_ns"]),
            rho=float(mat["rho_kg_mm3"]),
            cp=float(mat["cp_J_kgK"]),
            k=float(mat["k_kW_mmK"]),
            T_init=float(mat["T_init_K"]),
        )
        k_name = f"v17_{case.name}.k"
        k_path = args.out_dir / k_name
        nbytes = write_k(build_k_text(case, mesh, cfg), k_path)
        print(f"[OK] {k_path.name}  ({nbytes:,} bytes, Ep={case.Ep_uJ} uJ, "
              f"I0={case.I0_peak_kW_per_mm2:.3e} kW/mm^2)")

        registry_rows.append({
            "name": case.name,
            "Ep_uJ": f"{case.Ep_uJ:g}",
            "tau_ps": f"{case.tau_ps:g}",
            "w0_um": f"{case.w0_um:g}",
            "A": f"{case.A:g}",
            "I0_peak_kW_mm2": f"{case.I0_peak_kW_per_mm2:.6e}",
            "t_end_ns": f"{case.t_end_ns:g}",
            "k_file": f"models/v17_30ps_local/{k_name}",
            "run_lsdyna_final": "true" if case.run_lsdyna else "false",
            "NR": st["NR"], "NZ": st["NZ"],
            "n_nodes": st["n_nodes"], "n_elements": st["n_elements"],
            "dz_min_nm": f"{st['dz_min_nm']:.4f}",
            "axis_top_node_id": axis_top_node_id(mesh),
            "notes": case.notes,
        })

    reg_path = args.out_dir / "v17_case_registry.csv"
    with reg_path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(registry_rows[0].keys()))
        wr.writeheader()
        wr.writerows(registry_rows)
    print(f"\n[OK] registry -> {reg_path}")
    print(f"[DONE] {len(registry_rows)} cases written to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
