"""
plot_v50_ejecta_plume_evolution.py
==================================

Create paper-style V5A evolution figures directly from the real d3plot:
ejecta marker cloud, silicon displacement, and pressure/stress wave.

Run after the V5A LS-DYNA case finishes:
    python scripts\\plot_v50_ejecta_plume_evolution.py --case ep_5p0uj_ejecta_8500ms
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import Normalize, TwoSlopeNorm  # noqa: E402


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


def import_lasso():
    try:
        from lasso.dyna import ArrayType, D3plot  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on local install
        raise SystemExit(
            "Missing dependency: lasso-python\n"
            "Install it once with:\n"
            "    python -m pip install lasso-python\n"
        ) from exc
    return D3plot, ArrayType


def configure_fonts() -> None:
    candidates = [
        "Microsoft YaHei",
        "SimHei",
        "SimSun",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in candidates:
        if name in available:
            plt.rcParams["font.sans-serif"] = [name]
            break
    plt.rcParams["axes.unicode_minus"] = False


def get_array(arrays: dict, array_type, name: str):
    key = getattr(array_type, name, name)
    if key in arrays:
        return arrays[key]
    if name in arrays:
        return arrays[name]
    return None


def state_time_ns(state: int, dt_output_ns: float) -> float:
    return max(state - 1, 0) * dt_output_ns


def format_time_label(state: int, dt_output_ns: float) -> str:
    t_ns = state_time_ns(state, dt_output_ns)
    if t_ns < 1.0:
        return f"{t_ns * 1000.0:.0f} ps"
    return f"{t_ns:.2f} ns"


def pressure_from_shell_stress(shell_stress: np.ndarray) -> np.ndarray:
    return -(shell_stress[..., 0] + shell_stress[..., 1] + shell_stress[..., 2]) / 3.0


def shell_centroids(coords: np.ndarray, shell_node_indexes: np.ndarray) -> np.ndarray:
    idx = np.asarray(shell_node_indexes, dtype=np.int64)
    if idx.ndim == 1:
        idx = idx.reshape((-1, 4))
    return coords[idx[:, :4], :2].mean(axis=1)


def read_ejecta_node_ids(k_path: Path) -> np.ndarray:
    lines = k_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    title_idx = None
    for i, line in enumerate(lines):
        if "Ejecta marker plume" in line:
            title_idx = i
            break
    if title_idx is None:
        raise SystemExit(f"could not find ejecta marker part in {k_path}")
    node_start = None
    for i in range(title_idx, len(lines)):
        if lines[i].strip().upper() == "*NODE":
            node_start = i + 1
            break
    if node_start is None:
        raise SystemExit(f"could not find ejecta *NODE block in {k_path}")

    ids: list[int] = []
    for line in lines[node_start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("$"):
            continue
        if stripped.startswith("*"):
            break
        try:
            ids.append(int(stripped.split()[0]))
        except Exception:
            pass
    if not ids:
        raise SystemExit(f"no ejecta node ids parsed from {k_path}")
    return np.asarray(ids, dtype=np.int64)


def mirror_axisymmetric(x_um: np.ndarray, y_um: np.ndarray,
                        *values: np.ndarray) -> tuple[np.ndarray, ...]:
    keep = x_um > 1e-12
    out = [np.concatenate([-x_um[keep], x_um]), np.concatenate([y_um[keep], y_um])]
    for value in values:
        out.append(np.concatenate([value[keep], value]))
    return tuple(out)


def crop_mask(x_um: np.ndarray, y_um: np.ndarray, xlim: tuple[float, float],
              ylim: tuple[float, float]) -> np.ndarray:
    return (x_um >= xlim[0]) & (x_um <= xlim[1]) & (y_um >= ylim[0]) & (y_um <= ylim[1])


def levels(norm: Normalize, count: int = 18) -> np.ndarray:
    return np.linspace(float(norm.vmin), float(norm.vmax), count + 1)


def add_direction_marker(ax) -> None:
    x0, y0 = -0.06, -0.42
    dx, dy = 0.16, 0.24
    trans = ax.transAxes
    ax.annotate("", xy=(x0 + dx, y0), xytext=(x0, y0),
                xycoords=trans, textcoords=trans, annotation_clip=False,
                arrowprops={"arrowstyle": "->", "lw": 1.4, "color": "black"})
    ax.annotate("", xy=(x0, y0 + dy), xytext=(x0, y0),
                xycoords=trans, textcoords=trans, annotation_clip=False,
                arrowprops={"arrowstyle": "->", "lw": 1.4, "color": "black"})
    ax.text(x0 + dx * 1.08, y0, "r", va="center", ha="left", fontsize=8,
            transform=trans, clip_on=False)
    ax.text(x0, y0 + dy * 1.10, "z", va="bottom", ha="center", fontsize=8,
            transform=trans, clip_on=False)
    ax.text(x0 + dx * 0.12, y0 - 0.10, "Si section", va="top", ha="left",
            fontsize=7, transform=trans, clip_on=False)


def plot_panels(case: str, states: list[int], dt_output_ns: float, out_dir: Path,
                coords: np.ndarray, node_ids: np.ndarray, current: np.ndarray,
                velocity: np.ndarray | None, shell_idx: np.ndarray,
                shell_stress: np.ndarray, ejecta_node_ids: np.ndarray,
                xlim_um: tuple[float, float], zlim_um: tuple[float, float],
                plume_zlim_um: tuple[float, float], deform_scale: float) -> tuple[Path, Path]:
    ejecta_mask = np.isin(node_ids, ejecta_node_ids)
    substrate_node_mask = ~ejecta_mask

    cent = shell_centroids(coords, shell_idx) * 1000.0
    substrate_shell_mask = cent[:, 1] <= zlim_um[1] + 0.25

    disp = current - coords[None, :, :]
    disp_nm = np.linalg.norm(disp[:, substrate_node_mask, :2], axis=2) * 1.0e6
    disp_norm = Normalize(vmin=0.0, vmax=max(float(np.nanpercentile(disp_nm, 99.2)), 1.0))

    pressure = pressure_from_shell_stress(shell_stress)
    if pressure.ndim == 4:
        pressure = pressure.mean(axis=2)
    elif pressure.ndim == 3:
        pressure = pressure.mean(axis=2)
    p_sub = pressure[:, substrate_shell_mask]
    p_abs = max(float(np.nanpercentile(np.abs(p_sub), 99.0)), 1.0e-6)
    p_norm = TwoSlopeNorm(vmin=-p_abs, vcenter=0.0, vmax=p_abs)

    if velocity is not None:
        ejecta_speed = np.linalg.norm(velocity[:, ejecta_mask, :2], axis=2)
        plume_color_norm = Normalize(vmin=0.0, vmax=max(float(np.nanpercentile(ejecta_speed, 98.0)), 1.0))
    else:
        ejecta_speed = None
        plume_color_norm = Normalize(vmin=0.0, vmax=plume_zlim_um[1])

    cols = len(states)
    fig, axes = plt.subplots(3, cols, figsize=(3.65 * cols, 8.2),
                             dpi=220, sharex=False, constrained_layout=True)

    plume_sc = None
    disp_cn = None
    p_cn = None
    maxima_rows = []

    for j, state in enumerate(states):
        local_i = j
        label = format_time_label(state, dt_output_ns)

        cur_um = current[local_i] * 1000.0
        base_um = coords * 1000.0
        d_um = (current[local_i] - coords) * 1000.0

        # Row 1: ejecta/plume marker cloud.
        ax = axes[0, j]
        ex = cur_um[ejecta_mask, 0]
        ez = cur_um[ejecta_mask, 1]
        if velocity is not None:
            c = ejecta_speed[local_i]
            plume_sc = ax.scatter(ex, ez, c=c, s=8, cmap="hot", norm=plume_color_norm, linewidths=0)
        else:
            plume_sc = ax.scatter(ex, ez, c=ez, s=8, cmap="hot", norm=plume_color_norm, linewidths=0)
        ax.axvline(0.0, color="0.55", lw=0.8)
        ax.axhline(5.0, color="0.35", lw=0.9)
        ax.set_xlim(*xlim_um)
        ax.set_ylim(*plume_zlim_um)
        ax.set_title(label, fontsize=9)
        ax.set_aspect("equal", adjustable="box")

        # Row 2: deformed silicon displacement, mirrored around symmetry axis.
        ax = axes[1, j]
        sx = base_um[substrate_node_mask, 0]
        sz = base_um[substrate_node_mask, 1]
        sux = d_um[substrate_node_mask, 0]
        suz = d_um[substrate_node_mask, 1]
        smag = np.linalg.norm(d_um[substrate_node_mask, :2], axis=1) * 1000.0
        sx, sz, sux, suz, smag = mirror_axisymmetric(sx, sz, sux, suz, smag)
        mask = crop_mask(sx, sz, (-xlim_um[1], xlim_um[1]), zlim_um)
        tri = mtri.Triangulation(sx[mask] + sux[mask] * deform_scale,
                                 sz[mask] + suz[mask] * deform_scale)
        disp_cn = ax.tricontourf(tri, smag[mask], levels=levels(disp_norm), cmap="turbo", norm=disp_norm)
        ax.axvline(0.0, color="white", lw=0.8, alpha=0.8)
        ax.set_xlim(-xlim_um[1], xlim_um[1])
        ax.set_ylim(*zlim_um)
        ax.set_aspect("equal", adjustable="box")

        # Row 3: pressure/stress wave in silicon, mirrored around symmetry axis.
        ax = axes[2, j]
        px = cent[substrate_shell_mask, 0]
        pz = cent[substrate_shell_mask, 1]
        pp = p_sub[local_i]
        px, pz, pp = mirror_axisymmetric(px, pz, pp)
        mask = crop_mask(px, pz, (-xlim_um[1], xlim_um[1]), zlim_um)
        tri = mtri.Triangulation(px[mask], pz[mask])
        p_cn = ax.tricontourf(tri, pp[mask], levels=levels(p_norm), cmap="turbo", norm=p_norm)
        ax.axvline(0.0, color="white", lw=0.8, alpha=0.8)
        ax.set_xlim(-xlim_um[1], xlim_um[1])
        ax.set_ylim(*zlim_um)
        ax.set_aspect("equal", adjustable="box")

        maxima_rows.append((
            state,
            state_time_ns(state, dt_output_ns),
            float(np.nanmax(ez - 5.0)),
            float(np.nanmax(smag)),
            float(np.nanmax(pp)),
            float(np.nanmin(pp)),
        ))

    axes[0, 0].set_ylabel("ejecta z / um")
    axes[1, 0].set_ylabel("Si z / um")
    axes[2, 0].set_ylabel("Si z / um")
    for ax in axes[2, :]:
        ax.set_xlabel("signed r / um")
    for ax in axes[0, :]:
        ax.set_xlabel("r / um")
    add_direction_marker(axes[2, 0])

    fig.suptitle("5 µJ 激光烧蚀后硅片形变、应力波与喷射物演化图（0-5 ns）", fontsize=15)
    if plume_sc is not None:
        label = "ejecta speed / mm ms$^{-1}$" if velocity is not None else "ejecta height / um"
        fig.colorbar(plume_sc, ax=axes[0, :].tolist(), shrink=0.82, label=label)
    if disp_cn is not None:
        fig.colorbar(disp_cn, ax=axes[1, :].tolist(), shrink=0.82, label="|u| / nm")
    if p_cn is not None:
        fig.colorbar(p_cn, ax=axes[2, :].tolist(), shrink=0.82, label="pressure / GPa")

    out_dir.mkdir(parents=True, exist_ok=True)
    panel_path = out_dir / "v50_ejecta_plume_evolution_panels.png"
    fig.savefig(panel_path, bbox_inches="tight")
    plt.close(fig)

    arr = np.asarray(maxima_rows, dtype=float)
    fig2, ax1 = plt.subplots(figsize=(8.4, 4.8), dpi=180)
    ax2 = ax1.twinx()
    ax1.plot(arr[:, 1], arr[:, 2], "o-", color="#d62728", label="max ejecta height")
    ax1.plot(arr[:, 1], arr[:, 3], "s-", color="#1f77b4", label="max Si displacement")
    ax2.plot(arr[:, 1], arr[:, 4], "^-", color="#ff7f0e", label="max pressure")
    ax2.plot(arr[:, 1], arr[:, 5], "v-", color="#2ca02c", label="min pressure")
    ax1.set_xlabel("time / ns")
    ax1.set_ylabel("height or displacement / um, nm")
    ax2.set_ylabel("pressure / GPa")
    ax1.grid(True, alpha=0.3)
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="best")
    ax1.set_title("V5A real d3plot maxima over time")
    curve_path = out_dir / "v50_ejecta_plume_maxima.png"
    fig2.savefig(curve_path, bbox_inches="tight")
    plt.close(fig2)

    csv_path = out_dir / "v50_ejecta_plume_summary.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("state,time_ns,max_ejecta_height_above_surface_um,max_si_displacement_nm,max_pressure_GPa,min_pressure_GPa\n")
        for row in maxima_rows:
            fh.write(",".join(f"{v:.9g}" if isinstance(v, float) else str(v) for v in row) + "\n")
    return panel_path, curve_path


def main() -> int:
    configure_fonts()
    default_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=default_root)
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="ep_5p0uj_ejecta_8500ms")
    ap.add_argument("--states", default="1,50,100,200,377,502")
    ap.add_argument("--dt-output-ns", type=float, default=0.01)
    ap.add_argument("--xmax-um", type=float, default=25.0)
    ap.add_argument("--zmax-um", type=float, default=5.0)
    ap.add_argument("--plume-zmax-um", type=float, default=55.0)
    ap.add_argument("--deform-scale", type=float, default=200.0)
    args = ap.parse_args()

    root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else root
    case = args.case
    states = [int(v.strip()) for v in args.states.split(",") if v.strip()]

    d3plot_path = root / "results" / "v50_dynamic_ejecta_particle_pilot" / case / "d3plot"
    k_path = root / "models" / "v50_dynamic_ejecta_particle_pilot" / f"v50_{case}.k"
    if not d3plot_path.is_file():
        raise SystemExit(f"d3plot not found: {d3plot_path}")
    if not k_path.is_file():
        raise SystemExit(f"V50 keyword file not found: {k_path}")

    D3plot, ArrayType = import_lasso()
    d3 = D3plot(
        str(d3plot_path),
        state_filter={s - 1 for s in states},
        state_array_filter=[
            "node_displacement",
            "node_velocity",
            "element_shell_stress",
        ],
    )
    arrays = d3.arrays
    node_ids = get_array(arrays, ArrayType, "node_ids")
    coords = get_array(arrays, ArrayType, "node_coordinates")
    current = get_array(arrays, ArrayType, "node_displacement")
    velocity = get_array(arrays, ArrayType, "node_velocity")
    shell_idx = get_array(arrays, ArrayType, "element_shell_node_indexes")
    shell_stress = get_array(arrays, ArrayType, "element_shell_stress")
    if node_ids is None or coords is None or current is None or shell_idx is None or shell_stress is None:
        raise SystemExit("d3plot did not contain the required V50 arrays")

    ejecta_node_ids = read_ejecta_node_ids(k_path)
    out_dir = output_root / "figures" / "v50_ejecta_plume_evolution" / case
    panel, curve = plot_panels(
        case=case,
        states=states,
        dt_output_ns=args.dt_output_ns,
        out_dir=out_dir,
        coords=coords,
        node_ids=node_ids,
        current=current,
        velocity=velocity,
        shell_idx=shell_idx,
        shell_stress=shell_stress,
        ejecta_node_ids=ejecta_node_ids,
        xlim_um=(0.0, args.xmax_um),
        zlim_um=(0.0, args.zmax_um),
        plume_zlim_um=(args.zmax_um, args.plume_zmax_um),
        deform_scale=args.deform_scale,
    )
    print(f"[OK] {panel}")
    print(f"[OK] {curve}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
