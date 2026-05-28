"""
plot_v45_real_field_evolution.py
================================

Create paper-style multi-frame evolution figures from V44 real field CSV files.

Run after:
    python scripts\\extract_v44_lasso_d3plot.py --case ep_5p0uj_recoil_5000ms

Then:
    python scripts\\plot_v45_real_field_evolution.py --case ep_5p0uj_recoil_5000ms
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import Normalize, TwoSlopeNorm  # noqa: E402
from matplotlib import font_manager  # noqa: E402


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


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


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def load_node_field(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = read_csv(path)
    x = np.array([float(r["x_mm"]) * 1000.0 for r in rows])
    y = np.array([float(r["y_mm"]) * 1000.0 for r in rows])
    ux_raw = np.array([float(r["ux_mm"]) * 1000.0 for r in rows])
    uy_raw = np.array([float(r["uy_mm"]) * 1000.0 for r in rows])
    # Older V44 CSVs stored current coordinates in ux/uy.  Corrected V44 CSVs
    # store true displacement.  Detect and support both, so old outputs can
    # still be replotted without re-reading the d3plot.
    if np.nanmax(np.abs(ux_raw)) > 0.25 * max(np.nanmax(np.abs(x)), 1e-12):
        ux = ux_raw - x
        uy = uy_raw - y
    else:
        ux = ux_raw
        uy = uy_raw
    mag_nm = np.sqrt(ux * ux + uy * uy) * 1000.0
    return x, y, ux, uy, mag_nm


def load_pressure(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = read_csv(path)
    x = np.array([float(r["x_mm"]) * 1000.0 for r in rows])
    y = np.array([float(r["y_mm"]) * 1000.0 for r in rows])
    p = np.array([float(r["pressure_GPa"]) for r in rows])
    return x, y, p


def state_dirs(field_root: Path, case: str) -> list[Path]:
    case_dir = field_root / case
    return sorted([p for p in case_dir.glob("state_*") if p.is_dir()])


def pick_states(dirs: list[Path], wanted: list[int]) -> list[Path]:
    by_state = {int(p.name.split("_")[-1]): p for p in dirs}
    return [by_state[s] for s in wanted if s in by_state]


def has_field_csvs(state_dir: Path) -> bool:
    return (state_dir / "shell_pressure.csv").is_file() and (state_dir / "node_displacement.csv").is_file()


def state_number(state_dir: Path) -> int:
    return int(state_dir.name.split("_")[-1])


def state_time_ns(state: int, dt_output_ns: float = 0.01) -> float:
    # LS-DYNA writes state 1 at t=0, then state 2 at one output interval.
    return max(state - 1, 0) * dt_output_ns


def format_time_label(state: int, dt_output_ns: float = 0.01) -> str:
    t_ns = state_time_ns(state, dt_output_ns)
    if t_ns < 1.0:
        return f"{t_ns * 1000.0:.0f} ps"
    return f"{t_ns:.2f} ns"


def mirror_axisymmetric(x_um: np.ndarray, y_um: np.ndarray,
                        *values: np.ndarray) -> tuple[np.ndarray, ...]:
    keep = x_um > 1e-12
    mirrored = [np.concatenate([-x_um[keep], x_um])]
    mirrored.append(np.concatenate([y_um[keep], y_um]))
    for value in values:
        mirrored.append(np.concatenate([value[keep], value]))
    return tuple(mirrored)


def crop_mask(x_um: np.ndarray, y_um: np.ndarray, xlim: tuple[float, float],
              ylim: tuple[float, float]) -> np.ndarray:
    return (x_um >= xlim[0]) & (x_um <= xlim[1]) & (y_um >= ylim[0]) & (y_um <= ylim[1])


def fixed_contour_levels(norm: Normalize, count: int = 17) -> np.ndarray:
    """Return global contour boundaries so panels do not autoscale independently."""
    return np.linspace(float(norm.vmin), float(norm.vmax), count + 1)


def draw_pressure_panel(ax, state_dir: Path, xlim: tuple[float, float], ylim: tuple[float, float],
                        norm: TwoSlopeNorm, centered: bool, dt_output_ns: float):
    x, y, p = load_pressure(state_dir / "shell_pressure.csv")
    if centered:
        x, y, p = mirror_axisymmetric(x, y, p)
    mask = crop_mask(x, y, xlim, ylim)
    tri = mtri.Triangulation(x[mask], y[mask])
    cn = ax.tricontourf(tri, p[mask], levels=fixed_contour_levels(norm), cmap="turbo", norm=norm)
    ax.axvline(0.0, color="white", lw=0.8, alpha=0.75)
    ax.plot([0.0], [ylim[1]], marker="v", ms=3.5, color="black", clip_on=False)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(format_time_label(state_number(state_dir), dt_output_ns), fontsize=9)
    return cn


def draw_deformed_panel(ax, state_dir: Path, xlim: tuple[float, float], ylim: tuple[float, float],
                        magnification: float, centered: bool, dt_output_ns: float, norm=None):
    x, y, ux, uy, mag_nm = load_node_field(state_dir / "node_displacement.csv")
    if centered:
        x, y, ux, uy, mag_nm = mirror_axisymmetric(x, y, ux, uy, mag_nm)
    mask = crop_mask(x, y, xlim, ylim)
    xd = x[mask] + ux[mask] * magnification
    yd = y[mask] + uy[mask] * magnification
    values = mag_nm[mask]
    tri = mtri.Triangulation(xd, yd)
    cn = ax.tricontourf(tri, values, levels=fixed_contour_levels(norm), cmap="turbo", norm=norm)
    # Surface trace: original top row plus exaggerated deformed top row.
    top = mask & (np.abs(y - y.max()) < 1e-6)
    if top.any():
        order = np.argsort(x[top])
        ax.plot(x[top][order], y[top][order], color="white", lw=0.9, alpha=0.9)
        ax.plot((x[top] + ux[top] * magnification)[order],
                (y[top] + uy[top] * magnification)[order],
                color="black", lw=1.0)
    ax.axvline(0.0, color="white", lw=0.8, alpha=0.75)
    ax.plot([0.0], [ylim[1]], marker="v", ms=3.5, color="black", clip_on=False)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(format_time_label(state_number(state_dir), dt_output_ns), fontsize=9)
    return cn


def add_direction_marker(ax, xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
    x0 = -0.06
    y0 = -0.42
    dx = 0.16
    dy = 0.24
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


def as_axes_grid(axes) -> np.ndarray:
    axes_arr = np.asarray(axes, dtype=object)
    if axes_arr.ndim == 0:
        axes_arr = axes_arr.reshape((1, 1))
    elif axes_arr.ndim == 1:
        axes_arr = axes_arr.reshape((1, -1))
    return axes_arr


def hide_extra_axes(axes: np.ndarray, used: int) -> None:
    flat = axes.ravel()
    for ax in flat[used:]:
        ax.set_visible(False)


def plot_pressure_grid(dirs: list[Path], figure_root: Path, case: str,
                       xlim: tuple[float, float], ylim: tuple[float, float],
                       norm: TwoSlopeNorm, centered: bool, dt_output_ns: float,
                       rows: int) -> Path:
    cols = math.ceil(len(dirs) / rows)
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 2.45 * rows),
                             dpi=220, sharex=True, sharey=True,
                             constrained_layout=True)
    axes = as_axes_grid(axes)
    cn = None
    for ax, d in zip(axes.ravel(), dirs):
        cn = draw_pressure_panel(ax, d, xlim, ylim, norm, centered, dt_output_ns)
    hide_extra_axes(axes, len(dirs))
    add_direction_marker(axes[-1, 0], xlim, ylim)
    for ax in axes[:, 0]:
        ax.set_ylabel("z / um")
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel("signed r / um" if centered else "r / um")
    fig.suptitle("5 µJ 激光烧蚀后硅片内部反冲压力/应力波演化图（0–5 ns）", fontsize=14)
    if cn is not None:
        fig.colorbar(cn, ax=axes.ravel().tolist(), shrink=0.88, label="pressure / GPa")
    out = figure_root / "v45_pressure_evolution_grid.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_deformation_grid(dirs: list[Path], figure_root: Path, case: str,
                          xlim: tuple[float, float], ylim: tuple[float, float],
                          magnification: float, centered: bool, dt_output_ns: float,
                          rows: int, norm: Normalize) -> Path:
    cols = math.ceil(len(dirs) / rows)
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols, 2.45 * rows),
                             dpi=220, sharex=True, sharey=True,
                             constrained_layout=True)
    axes = as_axes_grid(axes)
    cn = None
    for ax, d in zip(axes.ravel(), dirs):
        cn = draw_deformed_panel(ax, d, xlim, ylim, magnification, centered, dt_output_ns, norm)
    hide_extra_axes(axes, len(dirs))
    add_direction_marker(axes[-1, 0], xlim, ylim)
    for ax in axes[:, 0]:
        ax.set_ylabel("z / um")
    for ax in axes[-1, :]:
        if ax.get_visible():
            ax.set_xlabel("signed r / um" if centered else "r / um")
    fig.suptitle("5 µJ 激光烧蚀后硅片在反冲压力作用下的位移演化图（0–5 ns）", fontsize=14)
    if cn is not None:
        fig.colorbar(cn, ax=axes.ravel().tolist(), shrink=0.88, label="|u| / nm")
    out = figure_root / "v45_deformation_evolution_grid.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def summary_curve(summary_csv: Path, case: str):
    rows = [r for r in read_csv(summary_csv) if r["case"] == case]
    by_field: dict[str, list[tuple[int, float]]] = {}
    for r in rows:
        by_field.setdefault(r["field"], []).append((int(r["state"]), float(r["max_value"])))
    return by_field


def main() -> int:
    configure_fonts()
    default_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=default_root)
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="ep_5p0uj_recoil_5000ms")
    ap.add_argument("--states", default="1,50,100,150,200,255,318,377,450,502")
    ap.add_argument("--dt-output-ns", type=float, default=0.01,
                    help="d3plot output interval; V42 default is 0.01 ns")
    ap.add_argument("--xmax-um", type=float, default=25.0)
    ap.add_argument("--ymax-um", type=float, default=5.0)
    ap.add_argument("--deform-scale", type=float, default=200.0)
    ap.add_argument("--grid-rows", type=int, default=2,
                    help="rows for separate pressure/deformation grid figures")
    ap.add_argument("--one-sided", action="store_true",
                    help="show the raw r>=0 section instead of mirrored centered view")
    args = ap.parse_args()

    root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else root
    field_root = root / "field_data" / "v44_lsprepost_real_field_extract"
    figure_root = output_root / "figures" / "v45_real_field_evolution" / args.case
    figure_root.mkdir(parents=True, exist_ok=True)
    states = [int(v.strip()) for v in args.states.split(",") if v.strip()]
    dirs = pick_states(state_dirs(field_root, args.case), states)
    if not dirs:
        raise SystemExit(f"no V44 field data found for {args.case}")
    missing = [d.name for d in dirs if not has_field_csvs(d)]
    dirs = [d for d in dirs if has_field_csvs(d)]
    if missing:
        print(f"[WARN] skipped states without field CSVs: {', '.join(missing)}")
    if not dirs:
        raise SystemExit(f"no usable V44 field CSVs found for {args.case}; run the V44 extractor first")

    centered = not args.one_sided
    xlim = (-args.xmax_um, args.xmax_um) if centered else (0.0, args.xmax_um)
    ylim = (0.0, args.ymax_um)
    n = len(dirs)

    all_p = []
    all_u = []
    for d in dirs:
        _, _, p = load_pressure(d / "shell_pressure.csv")
        all_p.append(p)
        *_, mag = load_node_field(d / "node_displacement.csv")
        all_u.append(mag)
    p_abs = max(float(max(np.nanmax(np.abs(v)) for v in all_p)), 1e-12)
    u_max = float(max(np.nanmax(v) for v in all_u))
    u_norm = Normalize(vmin=0.0, vmax=max(u_max, 1e-12))
    p_norm = TwoSlopeNorm(vmin=-p_abs, vcenter=0.0, vmax=p_abs)

    fig, axes = plt.subplots(2, n, figsize=(1.95 * n, 4.2), dpi=220, sharex=True, sharey=True)
    for ax, d in zip(axes[0], dirs):
        cn_p = draw_pressure_panel(ax, d, xlim, ylim, p_norm, centered, args.dt_output_ns)
    for ax, d in zip(axes[1], dirs):
        cn_u = draw_deformed_panel(ax, d, xlim, ylim, args.deform_scale, centered, args.dt_output_ns, u_norm)

    axes[0, 0].set_ylabel("pressure\nz / um")
    axes[1, 0].set_ylabel("deformed\nz / um")
    for ax in axes[-1]:
        ax.set_xlabel("signed r / um" if centered else "r / um")
    fig.colorbar(cn_p, ax=axes[0, :], shrink=0.78, label="pressure / GPa")
    fig.colorbar(cn_u, ax=axes[1, :], shrink=0.78, label="|u| / nm")
    fig.suptitle(f"{args.case}: real d3plot field evolution (axisymmetric section)", y=0.99)
    fig.savefig(figure_root / "v45_real_field_evolution_panels.png", bbox_inches="tight")
    if centered:
        fig.savefig(figure_root / "v45_centered_evolution_panels.png", bbox_inches="tight")
    plt.close(fig)

    pressure_grid = plot_pressure_grid(
        dirs, figure_root, args.case, xlim, ylim, p_norm,
        centered, args.dt_output_ns, max(1, args.grid_rows),
    )
    deformation_grid = plot_deformation_grid(
        dirs, figure_root, args.case, xlim, ylim, args.deform_scale,
        centered, args.dt_output_ns, max(1, args.grid_rows), u_norm,
    )

    curves = summary_curve(root / "lightweight_results" / "v44_lsprepost_real_field_extract" / "v44_real_field_summary.csv", args.case)
    fig, ax1 = plt.subplots(figsize=(7.2, 3.6), dpi=180)
    if "resultant_displacement" in curves:
        st, val = zip(*sorted(curves["resultant_displacement"]))
        ax1.plot(st, np.array(val) * 1.0e6, marker="o", label="max displacement")
        ax1.set_ylabel("max |u| / nm")
    ax2 = ax1.twinx()
    if "pressure_from_shell_stress" in curves:
        st, val = zip(*sorted(curves["pressure_from_shell_stress"]))
        ax2.plot(st, val, marker="s", color="tab:red", label="max pressure")
        ax2.set_ylabel("max pressure / GPa")
    ax1.set_xlabel("d3plot state")
    ax1.grid(True, alpha=0.25)
    fig.suptitle(f"{args.case}: real-field maxima over time")
    fig.savefig(figure_root / "v45_real_field_maxima.png", bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] {figure_root / 'v45_real_field_evolution_panels.png'}")
    print(f"[OK] {pressure_grid}")
    print(f"[OK] {deformation_grid}")
    print(f"[OK] {figure_root / 'v45_real_field_maxima.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
