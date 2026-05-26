"""
Build V1.7 grid-convergence LS-DYNA inputs for surface dz sensitivity.

This script generates focused comparison models for surface dz values such as
5 / 10 / 20 nm. It does not run LS-DYNA.

Run:
    python scripts\\build_v17_grid_convergence.py --dry-run
    python scripts\\build_v17_grid_convergence.py
"""

from __future__ import annotations

import argparse
import copy
import csv
import sys
import tomllib
from pathlib import Path

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from _build_v17_case import (  # noqa: E402
    V17Case,
    axis_top_node_id,
    build_k_text,
    build_mesh,
    resolve_database_dt_ms,
    write_k,
)
from build_v17_30ps_local_mesh import write_mesh_nodes_csv  # noqa: E402


def dz_label(dz_nm: float) -> str:
    text = f"{dz_nm:g}".replace(".", "p")
    return f"dz{text}nm"


def load_toml(path: Path) -> dict:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def select_cases(base_cfg: dict, names: set[str]) -> list[dict]:
    rows = [case for case in base_cfg["selected_cases"] if case["name"] in names]
    found = {case["name"] for case in rows}
    missing = sorted(names - found)
    if missing:
        raise ValueError(f"case(s) not found in base config: {', '.join(missing)}")
    return rows


def build_case(case_cfg: dict, mesh, cfg: dict) -> V17Case:
    mat = cfg["material"]
    laser = cfg["laser"]
    sol = cfg["solver"]
    return V17Case(
        name=str(case_cfg["name"]),
        Ep_uJ=float(case_cfg["Ep_uJ"]),
        run_lsdyna=bool(case_cfg["run_lsdyna_hint"]),
        notes=str(case_cfg.get("notes", "")),
        tau_ps=float(laser["pulse_width_ps"]),
        w0_um=float(laser["spot_radius_um"]),
        A=float(laser["absorption_A"]),
        t_end_ns=float(sol["t_end_ns"]),
        rho=float(mat["rho_kg_mm3"]),
        cp=float(mat["cp_J_kgK"]),
        k=float(mat["k_kW_mmK"]),
        T_init=float(mat["T_init_K"]),
    )


def write_registry(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--toml",
        type=Path,
        default=project_root / "config" / "v17_grid_convergence.toml",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    plan = load_toml(args.toml)
    base_cfg_path = project_root / plan["base"]["toml"]
    base_cfg = load_toml(base_cfg_path)
    out_root = project_root / plan["base"]["out_dir"]
    dz_values = [float(v) for v in plan["sweep"]["surface_dz_nm"]]
    case_names = set(str(v) for v in plan["sweep"]["case_names"])
    cases = select_cases(base_cfg, case_names)

    print("=" * 78)
    print("V1.7 grid convergence builder")
    print(f"  base config : {base_cfg_path}")
    print(f"  out root    : {out_root}")
    print(f"  dz values   : {', '.join(f'{v:g} nm' for v in dz_values)}")
    print(f"  cases       : {', '.join(case['name'] for case in cases)}")
    if args.dry_run:
        print("  mode        : dry-run")
    print("=" * 78)

    combined_rows: list[dict] = []
    for dz_nm in dz_values:
        cfg = copy.deepcopy(base_cfg)
        cfg["grid_z_um"]["surface_dz_nm"] = dz_nm
        cfg["selected_cases"] = cases
        mesh = build_mesh(cfg)
        st = mesh.stats
        label = dz_label(dz_nm)
        out_dir = out_root / label
        t_end_ms = float(cfg["solver"]["t_end_ns"]) * 1e-6
        dt_ms, dt_desc = resolve_database_dt_ms(cfg, t_end_ms)

        print()
        print(f"[{label}] NR x NZ = {st['NR']} x {st['NZ']}; "
              f"nodes={st['n_nodes']}; elements={st['n_elements']}; "
              f"dz_min={st['dz_min_nm']:.3f} nm; database DT={dt_ms * 1e6:g} ns")
        print(f"  {dt_desc}")

        if args.dry_run:
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        write_mesh_nodes_csv(mesh, out_dir / "mesh_nodes.csv")
        with (out_dir / "mesh_summary.csv").open("w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=[
                "surface_dz_nm", "NR", "NZ", "n_nodes", "n_elements",
                "dr_min_um", "dr_max_um", "dz_min_nm", "dz_max_nm",
                "axis_top_node_id",
            ])
            wr.writeheader()
            wr.writerow({
                "surface_dz_nm": f"{dz_nm:g}",
                "NR": st["NR"],
                "NZ": st["NZ"],
                "n_nodes": st["n_nodes"],
                "n_elements": st["n_elements"],
                "dr_min_um": f"{st['dr_min_um']:.6f}",
                "dr_max_um": f"{st['dr_max_um']:.6f}",
                "dz_min_nm": f"{st['dz_min_nm']:.4f}",
                "dz_max_nm": f"{st['dz_max_um'] * 1000:.4f}",
                "axis_top_node_id": axis_top_node_id(mesh),
            })

        for case_cfg in cases:
            case = build_case(case_cfg, mesh, cfg)
            k_name = f"v17_{label}_{case.name}.k"
            k_path = out_dir / k_name
            if not args.dry_run:
                nbytes = write_k(build_k_text(case, mesh, cfg), k_path)
                print(f"  [OK] {k_path.relative_to(project_root)} ({nbytes:,} bytes)")
            combined_rows.append({
                "grid_label": label,
                "surface_dz_nm": f"{dz_nm:g}",
                "case": case.name,
                "Ep_uJ": f"{case.Ep_uJ:g}",
                "I0_peak_kW_mm2": f"{case.I0_peak_kW_per_mm2:.6e}",
                "k_file": (out_dir / k_name).relative_to(project_root).as_posix(),
                "run_lsdyna_final": "true",
                "NR": st["NR"],
                "NZ": st["NZ"],
                "n_nodes": st["n_nodes"],
                "n_elements": st["n_elements"],
                "axis_top_node_id": axis_top_node_id(mesh),
                "notes": case.notes,
            })

    if not args.dry_run:
        write_registry(combined_rows, out_root / "grid_convergence_registry.csv")
        print()
        print(f"[DONE] registry -> {out_root / 'grid_convergence_registry.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
