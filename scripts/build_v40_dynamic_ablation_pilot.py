"""
build_v40_dynamic_ablation_pilot.py
===================================

Generate V4.0A structural-dynamics pilot LS-DYNA decks.

This is the first step toward dynamic ablation.  It reuses the V1.7 local mesh,
estimates the vapor-zone radius from the V1.7 vapor threshold, and assigns a
small ejecta impulse to near-surface nodes inside that radius.

Run from project root:
    python scripts\\build_v40_dynamic_ablation_pilot.py
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from _build_v17_case import (  # noqa: E402
    axis_top_node_id,
    build_mesh,
    f10,
    gen_elements,
    gen_nodes,
    i10,
    node_id_at,
    write_k,
)


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


@dataclass(frozen=True)
class V40Case:
    name: str
    Ep_uJ: float
    notes: str


def load_toml(path: Path) -> dict:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def vapor_radius_um(Ep_uJ: float, threshold_uJ: float, w0_um: float) -> float:
    if Ep_uJ <= threshold_uJ:
        return 0.0
    return w0_um / math.sqrt(2.0) * math.sqrt(math.log(Ep_uJ / threshold_uJ))


def top_layer_element_ids(mesh, radius_um: float, depth_um: float) -> list[int]:
    ids: list[int] = []
    if radius_um <= 0.0 or depth_um <= 0.0:
        return ids

    z_min = mesh.z_um[-1] - depth_um
    for j in range(mesh.NZ):
        zc = 0.5 * (mesh.z_um[j] + mesh.z_um[j + 1])
        if zc < z_min:
            continue
        for i in range(mesh.NR):
            rc = 0.5 * (mesh.r_um[i] + mesh.r_um[i + 1])
            if rc <= radius_um:
                ids.append(i + j * mesh.NR + 1)
    return ids


def ejecta_node_velocities(mesh, radius_um: float, depth_um: float,
                           vy_mm_ms: float, radial_fraction: float) -> list[str]:
    lines: list[str] = []
    if radius_um <= 0.0 or vy_mm_ms <= 0.0:
        return lines

    z_min = mesh.z_um[-1] - depth_um
    seen: set[int] = set()
    for j in range(mesh.NZ + 1):
        if mesh.z_um[j] < z_min:
            continue
        for i in range(mesh.NR + 1):
            r_um = float(mesh.r_um[i])
            if r_um > radius_um:
                continue
            nid = node_id_at(mesh, i, j)
            if nid in seen:
                continue
            seen.add(nid)
            profile = math.exp(-2.0 * (r_um / max(radius_um, 1e-9)) ** 2)
            vy = vy_mm_ms * profile
            vx = radial_fraction * vy * (r_um / max(radius_um, 1e-9))
            lines.append(
                f"{nid:8d}{vx:16.7e}{vy:16.7e}{0.0:16.7e}"
                f"{0.0:16.7e}{0.0:16.7e}{0.0:16.7e}"
            )
    return lines


def boundary_spc_lines(mesh) -> list[str]:
    lines: list[str] = []
    # Axis: fix radial displacement (x) and out-of-plane z for 2D axisymmetric use.
    for j in range(mesh.NZ + 1):
        nid = node_id_at(mesh, 0, j)
        lines.append(f"{nid:8d}{0:8d}{1:8d}{0:8d}{1:8d}{0:8d}{0:8d}{0:8d}")
    # Bottom: fix vertical displacement (y) and out-of-plane z.
    for i in range(mesh.NR + 1):
        nid = node_id_at(mesh, i, 0)
        lines.append(f"{nid:8d}{0:8d}{0:8d}{1:8d}{1:8d}{0:8d}{0:8d}{0:8d}")
    return lines


def build_dynamic_k_text(case: V40Case, mesh, cfg: dict, v17_cfg: dict) -> tuple[str, dict]:
    mat = cfg["material"]
    dyn = cfg["dynamic_solver"]
    abl = cfg["ablation_model"]
    laser = v17_cfg["laser"]

    t_end_ms = float(dyn["t_end_ns"]) * 1.0e-6
    dt_plot_ms = float(dyn["dt_plot_ns"]) * 1.0e-6
    threshold = float(abl["vapor_threshold_Ep_uJ"])
    w0_um = float(laser["spot_radius_um"])
    r_vap_um = vapor_radius_um(case.Ep_uJ, threshold, w0_um)
    depth_um = float(abl["active_layer_depth_um"])

    # Unit conversion: 1 m/s = 1 mm/ms.
    ejecta_velocity = float(abl["ejecta_velocity_m_s"])
    radial_fraction = float(abl["radial_velocity_fraction"])
    velocity_lines = ejecta_node_velocities(
        mesh=mesh,
        radius_um=r_vap_um,
        depth_um=depth_um,
        vy_mm_ms=ejecta_velocity,
        radial_fraction=radial_fraction,
    )
    erosion_eids = top_layer_element_ids(mesh, radius_um=r_vap_um, depth_um=depth_um)

    nodes = gen_nodes(mesh)
    elems = gen_elements(mesh, pid=1)
    spc = boundary_spc_lines(mesh)

    title = f"V40A dynamic ablation pilot {case.name} Ep={case.Ep_uJ:g}uJ"[:80]
    mat_line = (
        i10(1) + f10(float(mat["rho_kg_mm3"]))
        + f10(float(mat["youngs_modulus_GPa"]))
        + f10(float(mat["poisson_ratio"]))
    )
    sec_card = i10(1) + i10(15)
    sec_card2 = f10(1.0) + f10(1.0) + f10(1.0) + f10(1.0)
    part_card = i10(1) + i10(1) + i10(1)

    velocity_block = ""
    if velocity_lines:
        velocity_block = (
            "*INITIAL_VELOCITY_NODE\n"
            "$#     nid              vx              vy              vz"
            "             vrx             vry             vrz\n"
            + "\n".join(velocity_lines)
            + "\n$\n"
        )
    else:
        velocity_block = (
            "$ No initial ejecta velocity: case is at or below vapor threshold.\n$\n"
        )

    text = f"""$ ============================================================================
$  v40_{case.name}.k
$  Generated by scripts/build_v40_dynamic_ablation_pilot.py -- DO NOT HAND-EDIT.
$
$  V4.0A dynamic ablation pilot.
$  This is a structural dynamics scaffold, not a final calibrated ablation deck.
$
$  Source mesh: V1.7 local refined 30 ps mesh.
$  Case: {case.name}, Ep = {case.Ep_uJ:g} uJ
$  V1.7 vapor threshold estimate: {threshold:g} uJ
$  Estimated vapor radius: {r_vap_um:.4f} um
$  Active layer depth: {depth_um:g} um
$  Initial ejecta velocity: {ejecta_velocity:g} m/s
$
$  Units: mm - ms - kg - kN - GPa
$ ============================================================================
*KEYWORD
*TITLE
{title}
$
*CONTROL_TERMINATION
{f10(t_end_ms)}
$
*CONTROL_TIMESTEP
{f10(0.0)}{f10(float(dyn["timestep_scale"]))}
$
*DATABASE_BINARY_D3PLOT
{f10(dt_plot_ms)}
$
*DATABASE_GLSTAT
{f10(dt_plot_ms)}
$
*DATABASE_MATSUM
{f10(dt_plot_ms)}
$
*MAT_ELASTIC
$#     mid        ro         e        pr
{mat_line}
$
*SECTION_SHELL
{sec_card}
{sec_card2}
$
*PART
Silicon V40A dynamic pilot {case.name}
{part_card}
$
*NODE
"""
    text += "\n".join(nodes) + "\n$\n"
    text += "*ELEMENT_SHELL\n" + "\n".join(elems) + "\n$\n"
    text += "*BOUNDARY_SPC_NODE\n"
    text += "$#     nid     cid    dofx    dofy    dofz   dofrx   dofry   dofrz\n"
    text += "\n".join(spc) + "\n$\n"
    text += velocity_block
    text += f"""$ ---------------------------------------------------------------------------
$ Active material deletion is NOT enabled in V4.0A.
$ Candidate deletion elements are exported to v40_erosion_candidates.csv.
$ Next step V4.0B should enable the LS-DYNA erosion option supported by the
$ local solver/material version and compare crater growth against V2.6.
$ ---------------------------------------------------------------------------
*END
"""

    meta = {
        "vapor_radius_um": r_vap_um,
        "active_layer_depth_um": depth_um,
        "ejecta_velocity_m_s": ejecta_velocity,
        "ejecta_nodes": len(velocity_lines),
        "erosion_candidate_elements": len(erosion_eids),
        "axis_top_node_id": axis_top_node_id(mesh),
    }
    return text, meta


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--toml",
        type=Path,
        default=project_root / "config" / "v40_dynamic_ablation_pilot.toml",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=project_root / "models" / "v40_dynamic_ablation_pilot",
    )
    args = ap.parse_args()

    cfg = load_toml(args.toml)
    v17_cfg = load_toml(project_root / cfg["source"]["v17_config"])
    mesh = build_mesh(v17_cfg)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("V4.0A -- build dynamic ablation pilot decks")
    print(f"  source mesh: NR x NZ = {mesh.NR} x {mesh.NZ}")
    print(f"  nodes = {mesh.stats['n_nodes']}, elements = {mesh.stats['n_elements']}")
    print(f"  axis-top node id = {axis_top_node_id(mesh)}")
    print("=" * 78)

    registry_rows = []
    erosion_rows = []
    for sc in cfg["selected_cases"]:
        case = V40Case(
            name=str(sc["name"]),
            Ep_uJ=float(sc["Ep_uJ"]),
            notes=str(sc.get("notes", "")),
        )
        text, meta = build_dynamic_k_text(case, mesh, cfg, v17_cfg)
        k_name = f"v40_{case.name}.k"
        k_path = args.out_dir / k_name
        nbytes = write_k(text, k_path)

        eids = top_layer_element_ids(
            mesh,
            radius_um=float(meta["vapor_radius_um"]),
            depth_um=float(meta["active_layer_depth_um"]),
        )
        for eid in eids:
            j = (eid - 1) // mesh.NR
            i = (eid - 1) % mesh.NR
            erosion_rows.append({
                "case": case.name,
                "element_id": eid,
                "i": i,
                "j": j,
                "r_center_um": f"{0.5 * (mesh.r_um[i] + mesh.r_um[i + 1]):.6f}",
                "z_center_um": f"{0.5 * (mesh.z_um[j] + mesh.z_um[j + 1]):.6f}",
                "reason": "inside estimated vapor radius and active layer",
            })

        registry_rows.append({
            "name": case.name,
            "Ep_uJ": f"{case.Ep_uJ:g}",
            "k_file": f"models/v40_dynamic_ablation_pilot/{k_name}",
            "vapor_radius_um": f"{meta['vapor_radius_um']:.6f}",
            "active_layer_depth_um": f"{meta['active_layer_depth_um']:.6f}",
            "ejecta_velocity_m_s": f"{meta['ejecta_velocity_m_s']:.6f}",
            "ejecta_nodes": meta["ejecta_nodes"],
            "erosion_candidate_elements": meta["erosion_candidate_elements"],
            "axis_top_node_id": meta["axis_top_node_id"],
            "run_lsdyna_final": "true",
            "notes": case.notes,
        })
        print(
            f"[OK] {k_path.name} ({nbytes:,} bytes): "
            f"r_vap={meta['vapor_radius_um']:.3f} um, "
            f"ejecta_nodes={meta['ejecta_nodes']}, "
            f"candidate_elements={meta['erosion_candidate_elements']}"
        )

    reg_path = args.out_dir / "v40_case_registry.csv"
    with reg_path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(registry_rows[0].keys()))
        wr.writeheader()
        wr.writerows(registry_rows)
    print(f"[OK] registry -> {reg_path}")

    cand_path = args.out_dir / "v40_erosion_candidates.csv"
    fieldnames = [
        "case", "element_id", "i", "j", "r_center_um", "z_center_um", "reason",
    ]
    with cand_path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        wr.writerows(erosion_rows)
    print(f"[OK] erosion candidates -> {cand_path}")
    print("[DONE] V4.0A dynamic ablation pilot decks generated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
