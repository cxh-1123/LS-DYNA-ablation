"""
build_v41_dynamic_ablation_deletion.py
======================================

Generate V4B LS-DYNA decks with vapor-zone elements removed from the mesh.

V4B is intentionally robust: it represents material removal by omitting the
candidate vaporized elements from the generated structural mesh.  This avoids
solver-version-specific active erosion cards while still producing a true
crater geometry for dynamic response.

Run from project root:
    python scripts\\build_v41_dynamic_ablation_deletion.py
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
class V41Case:
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


def element_nodes(mesh, eid: int) -> tuple[int, int, int, int]:
    i = (eid - 1) % mesh.NR
    j = (eid - 1) // mesh.NR
    n1 = i + j * (mesh.NR + 1) + 1
    n2 = (i + 1) + j * (mesh.NR + 1) + 1
    n3 = (i + 1) + (j + 1) * (mesh.NR + 1) + 1
    n4 = i + (j + 1) * (mesh.NR + 1) + 1
    return n1, n2, n3, n4


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


def gen_elements_without(mesh, removed: set[int], pid: int = 1) -> list[str]:
    out: list[str] = []
    for j in range(mesh.NZ):
        for i in range(mesh.NR):
            eid = i + j * mesh.NR + 1
            if eid in removed:
                continue
            n1 = i + j * (mesh.NR + 1) + 1
            n2 = (i + 1) + j * (mesh.NR + 1) + 1
            n3 = (i + 1) + (j + 1) * (mesh.NR + 1) + 1
            n4 = i + (j + 1) * (mesh.NR + 1) + 1
            out.append(f"{eid:8d}{pid:8d}{n1:8d}{n2:8d}{n3:8d}{n4:8d}")
    return out


def boundary_spc_lines(mesh) -> list[str]:
    lines: list[str] = []
    for j in range(mesh.NZ + 1):
        nid = node_id_at(mesh, 0, j)
        lines.append(f"{nid:8d}{0:8d}{1:8d}{0:8d}{1:8d}{0:8d}{0:8d}{0:8d}")
    for i in range(mesh.NR + 1):
        nid = node_id_at(mesh, i, 0)
        lines.append(f"{nid:8d}{0:8d}{0:8d}{1:8d}{1:8d}{0:8d}{0:8d}{0:8d}")
    return lines


def exposed_boundary_nodes(mesh, removed: set[int]) -> set[int]:
    removed_nodes: set[int] = set()
    kept_nodes: set[int] = set()
    for eid in range(1, mesh.NR * mesh.NZ + 1):
        nodes = set(element_nodes(mesh, eid))
        if eid in removed:
            removed_nodes.update(nodes)
        else:
            kept_nodes.update(nodes)
    return removed_nodes.intersection(kept_nodes)


def recoil_velocity_lines(mesh, node_ids: set[int], radius_um: float,
                          vy_mm_ms: float, radial_fraction: float) -> list[str]:
    lines: list[str] = []
    if not node_ids or radius_um <= 0.0 or vy_mm_ms == 0.0:
        return lines
    for nid in sorted(node_ids):
        idx = nid - 1
        i = idx % (mesh.NR + 1)
        j = idx // (mesh.NR + 1)
        r_um = float(mesh.r_um[i])
        if r_um > radius_um + 1.0e-9:
            continue
        profile = math.exp(-2.0 * (r_um / max(radius_um, 1e-9)) ** 2)
        vy = vy_mm_ms * profile
        vx = radial_fraction * abs(vy) * (r_um / max(radius_um, 1e-9))
        lines.append(
            i10(nid) + f10(vx) + f10(vy) + f10(0.0)
            + f10(0.0) + f10(0.0) + f10(0.0)
        )
    return lines


def build_k_text(case: V41Case, mesh, cfg: dict, v17_cfg: dict) -> tuple[str, dict, list[int]]:
    mat = cfg["material"]
    dyn = cfg["dynamic_solver"]
    delete = cfg["deletion_model"]
    laser = v17_cfg["laser"]

    t_end_ms = float(dyn["t_end_ns"]) * 1.0e-6
    dt_plot_ms = float(dyn["dt_plot_ns"]) * 1.0e-6
    threshold = float(delete["vapor_threshold_Ep_uJ"])
    w0_um = float(laser["spot_radius_um"])
    r_vap_um = vapor_radius_um(case.Ep_uJ, threshold, w0_um)
    depth_um = float(delete["active_layer_depth_um"])
    removed = set(top_layer_element_ids(mesh, r_vap_um, depth_um))

    nodes = gen_nodes(mesh)
    elems = gen_elements_without(mesh, removed)
    spc = boundary_spc_lines(mesh)

    exposed = exposed_boundary_nodes(mesh, removed)
    velocity_lines = recoil_velocity_lines(
        mesh,
        exposed,
        radius_um=r_vap_um,
        vy_mm_ms=float(delete["recoil_velocity_m_s"]),
        radial_fraction=float(delete["radial_velocity_fraction"]),
    )

    title = f"V41 V4B predeleted ablation {case.name} Ep={case.Ep_uJ:g}uJ"[:80]
    mat_line = (
        i10(1) + f10(float(mat["rho_kg_mm3"]))
        + f10(float(mat["youngs_modulus_GPa"]))
        + f10(float(mat["poisson_ratio"]))
    )
    part_card = i10(1) + i10(1) + i10(1)

    velocity_block = ""
    if velocity_lines:
        velocity_block = (
            "*INITIAL_VELOCITY_NODE\n"
            "$#     nid        vx        vy        vz       vrx       vry       vrz\n"
            + "\n".join(velocity_lines)
            + "\n$\n"
        )
    else:
        velocity_block = "$ No recoil velocity nodes for this case.\n$\n"

    text = f"""$ ============================================================================
$  v41_{case.name}.k
$  Generated by scripts/build_v41_dynamic_ablation_deletion.py -- DO NOT HAND-EDIT.
$
$  V4B dynamic ablation with pre-deleted vapor-zone elements.
$  Removed elements: {len(removed)}
$  Estimated vapor radius: {r_vap_um:.4f} um
$  Removed layer depth: {depth_um:g} um
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
{i10(1)}{i10(15)}
{f10(1.0)}{f10(1.0)}{f10(1.0)}{f10(1.0)}
$
*PART
Silicon V4B predeleted ablation {case.name}
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
    text += """$ ---------------------------------------------------------------------------
$ V4B pre-deletes vapor-zone elements by omitting them from *ELEMENT_SHELL.
$ V4C should replace this with solver-version-specific active erosion.
$ ---------------------------------------------------------------------------
*END
"""

    meta = {
        "vapor_radius_um": r_vap_um,
        "removed_elements": len(removed),
        "kept_elements": len(elems),
        "recoil_nodes": len(velocity_lines),
        "axis_top_node_id": axis_top_node_id(mesh),
    }
    return text, meta, sorted(removed)


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--toml",
        type=Path,
        default=project_root / "config" / "v41_dynamic_ablation_deletion.toml",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=project_root / "models" / "v41_dynamic_ablation_deletion",
    )
    args = ap.parse_args()

    cfg = load_toml(args.toml)
    v17_cfg = load_toml(project_root / cfg["source"]["v17_config"])
    mesh = build_mesh(v17_cfg)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("V4B -- build pre-deleted dynamic ablation decks")
    print(f"  source mesh: NR x NZ = {mesh.NR} x {mesh.NZ}")
    print(f"  nodes = {mesh.stats['n_nodes']}, elements = {mesh.stats['n_elements']}")
    print("=" * 78)

    registry_rows = []
    removed_rows = []
    for sc in cfg["selected_cases"]:
        case = V41Case(
            name=str(sc["name"]),
            Ep_uJ=float(sc["Ep_uJ"]),
            notes=str(sc.get("notes", "")),
        )
        text, meta, removed = build_k_text(case, mesh, cfg, v17_cfg)
        k_name = f"v41_{case.name}.k"
        k_path = args.out_dir / k_name
        nbytes = write_k(text, k_path)

        for eid in removed:
            j = (eid - 1) // mesh.NR
            i = (eid - 1) % mesh.NR
            removed_rows.append({
                "case": case.name,
                "element_id": eid,
                "i": i,
                "j": j,
                "r_center_um": f"{0.5 * (mesh.r_um[i] + mesh.r_um[i + 1]):.6f}",
                "z_center_um": f"{0.5 * (mesh.z_um[j] + mesh.z_um[j + 1]):.6f}",
                "reason": "pre-deleted vapor-zone element",
            })

        registry_rows.append({
            "name": case.name,
            "Ep_uJ": f"{case.Ep_uJ:g}",
            "k_file": f"models/v41_dynamic_ablation_deletion/{k_name}",
            "vapor_radius_um": f"{meta['vapor_radius_um']:.6f}",
            "removed_elements": meta["removed_elements"],
            "kept_elements": meta["kept_elements"],
            "recoil_nodes": meta["recoil_nodes"],
            "axis_top_node_id": meta["axis_top_node_id"],
            "run_lsdyna_final": "true",
            "notes": case.notes,
        })
        print(
            f"[OK] {k_path.name} ({nbytes:,} bytes): "
            f"removed={meta['removed_elements']}, "
            f"kept={meta['kept_elements']}, recoil_nodes={meta['recoil_nodes']}"
        )

    reg_path = args.out_dir / "v41_case_registry.csv"
    with reg_path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(registry_rows[0].keys()))
        wr.writeheader()
        wr.writerows(registry_rows)
    print(f"[OK] registry -> {reg_path}")

    rem_path = args.out_dir / "v41_removed_elements.csv"
    fieldnames = [
        "case", "element_id", "i", "j", "r_center_um", "z_center_um", "reason",
    ]
    with rem_path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        wr.writerows(removed_rows)
    print(f"[OK] removed elements -> {rem_path}")
    print("[DONE] V4B decks generated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
