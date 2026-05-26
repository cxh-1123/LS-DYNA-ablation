"""
plot_v26_30ps_threshold_ablation.py
===================================

Generate V2.6 figures from V1.7 tprint + V2.6 threshold metrics CSVs.

Run from project root (after extract):
    .\\.venv\\Scripts\\python.exe scripts\\plot_v26_30ps_threshold_ablation.py
"""

from __future__ import annotations

import argparse
import csv
import sys
import tomllib
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from plot_v1_temperature import parse_tprint  # noqa: E402
from check_v17_outputs import tprint_to_grid  # noqa: E402
from extract_v26_30ps_threshold_ablation import (  # noqa: E402
    MAIN_CASE,
    DEFAULT_CASE_ORDER,
    V17MeshGrid,
    load_mesh_grid,
)

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "savefig.dpi": 300,
})

CASE_STYLE = {
    "ep_0p5uj": {"color": "#1f77b4", "marker": "o", "label": "ep_0p5uj"},
    "ep_1p0uj": {"color": "#2ca02c", "marker": "s", "label": "ep_1p0uj"},
    "ep_2p0uj": {"color": "#ff7f0e", "marker": "^", "label": "ep_2p0uj"},
    "ep_5p0uj": {"color": "#d62728", "marker": "D", "label": "ep_5p0uj (main)"},
}

V26_BANNER = (
    "V2.6 -- 30 ps threshold-equivalent ablation\n"
    "driver: V1.7 local mesh  |  main case: ep_5p0uj"
)

EVOLUTION_TARGETS_PS = [0, 10, 30, 60, 100, 300, 1000, 2000, 5000]
EARLY_TIME_PS_MAX = 120.0
EARLY_PLOT_CASES = ["ep_1p0uj", "ep_2p0uj", "ep_5p0uj"]
ZOOM_R_UM = (-50.0, 50.0)
ZOOM_Z_UM = (4.5, 5.0)
Z_TOP_UM = 5.0
ZOOM_Z_SURFACE_UM = (4.94, 5.005)  # ~65 nm window -- makes ~42 nm crater visible
DEPTH_PROFILE_MAX_UM = 0.08        # peak profile: depth below top surface


def pick_nearest_time_index(
    times_ns: np.ndarray,
    target_ns: float,
    *,
    case_name: str = "",
) -> tuple[int, float, float]:
    """Local copy to avoid importing matplotlib-heavy plot_v17."""
    times_ns = np.asarray(times_ns, dtype=float)
    t_min, t_max = float(times_ns[0]), float(times_ns[-1])
    label = f" ({case_name})" if case_name else ""
    if target_ns > t_max:
        idx = int(times_ns.size - 1)
        return idx, float(times_ns[idx]), abs(float(times_ns[idx]) - target_ns)
    if target_ns < t_min:
        return 0, t_min, abs(t_min - target_ns)
    if target_ns == 0.0:
        return 0, t_min, abs(t_min - target_ns)
    pool = np.arange(times_ns.size, dtype=int)
    if t_min == 0.0 and times_ns.size > 1:
        pool = pool[1:]
    idx = int(pool[int(np.argmin(np.abs(times_ns[pool] - target_ns)))])
    actual = float(times_ns[idx])
    return idx, actual, abs(actual - target_ns)


def mirror_field(r_um: np.ndarray, T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r_full = np.concatenate([-r_um[:0:-1], r_um])
    T_full = np.concatenate([T[:, :0:-1], T], axis=1)
    return r_full, T_full


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


def read_summary_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_tprint_grids(tprint_path: Path, mesh: V17MeshGrid) -> tuple[np.ndarray, list[np.ndarray]]:
    times_ms, T_all = parse_tprint(tprint_path)
    times_ns = np.asarray(times_ms, dtype=float) * 1e6
    grids = [
        tprint_to_grid(T_all[k], mesh.node_map, mesh.NR, mesh.NZ)
        for k in range(T_all.shape[0])
    ]
    return times_ns, grids


def add_threshold_contours(ax, r_full, z_full, T_full, T_melt, T_vap, T_local_max):
    handles = []
    try:
        ax.contour(r_full, z_full, T_full, levels=[T_melt], colors="#00e5ff", linewidths=1.2)
        handles.append(Line2D([0], [0], color="#00e5ff", lw=1.2, label=f"T_melt={T_melt:.0f} K"))
    except Exception:
        pass
    if T_local_max >= T_vap:
        try:
            ax.contour(r_full, z_full, T_full, levels=[T_vap], colors="#ff00aa", linewidths=1.6)
            handles.append(Line2D([0], [0], color="#ff00aa", lw=1.6, label=f"T_vap={T_vap:.0f} K"))
        except Exception:
            pass
    ax.axvline(0.0, color="white", linestyle="--", linewidth=0.7, alpha=0.8)
    return handles


def plot_case_comparison_table(summary_rows: list[dict[str, str]], out_path: Path) -> None:
    headers = [
        "case", "Tmax_peak\n(K)", "t_peak\n(ps)",
        "max melt d\n(um)", "max melt r\n(um)",
        "max vapor d\n(um)", "max vapor r\n(um)",
        "max crater d\n(um)", "max crater r\n(um)",
        "final regime",
    ]
    keys = [
        "case", "Tmax_peak_K", "t_peak_ps",
        "max_melt_depth_um", "max_melt_radius_um",
        "max_vapor_depth_um", "max_vapor_radius_um",
        "max_crater_depth_um", "max_crater_radius_um",
        "final_regime",
    ]
    cells, colors = [], []
    for row in summary_rows:
        if row.get("status") not in ("OK", ""):
            cells.append([row.get("case", ""), "missing"] + [""] * (len(keys) - 2))
            colors.append("#eeeeee")
            continue
        cells.append([row.get(k, "") for k in keys])
        regime = row.get("final_regime", "")
        if regime == "vapor/ablation candidate":
            colors.append("#ffe2e2")
        elif regime == "melt only":
            colors.append("#fff5d6")
        elif regime == "missing":
            colors.append("#eeeeee")
        else:
            colors.append("#eef3ff")

    fig, ax = plt.subplots(figsize=(14.0, max(3.2, 1.3 + 0.55 * len(cells))))
    ax.axis("off")
    tbl = ax.table(cellText=cells, colLabels=headers, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1.0, 1.75)
    for j in range(len(headers)):
        tbl[(0, j)].set_facecolor("#333333")
        tbl[(0, j)].set_text_props(color="white", weight="bold")
    for i, c in enumerate(colors, start=1):
        for j in range(len(headers)):
            tbl[(i, j)].set_facecolor(c)
    ax.set_title(V26_BANNER, pad=14, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_metric_vs_time(
    case_metrics: dict[str, dict[str, np.ndarray]],
    crater_key: str,
    melt_key: str,
    ylabel: str,
    out_path: Path,
    *,
    t_xlim_ps: tuple[float, float] | None = None,
    title_suffix: str = "",
    ep5_markers_ps: tuple[float, float] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    legend_crater, legend_melt = [], []
    for case in (EARLY_PLOT_CASES if t_xlim_ps else DEFAULT_CASE_ORDER):
        if case not in case_metrics:
            continue
        d = case_metrics[case]
        st = CASE_STYLE[case]
        t_ps = d["time_ns"] * 1e3
        if t_xlim_ps:
            mask = (t_ps >= t_xlim_ps[0]) & (t_ps <= t_xlim_ps[1])
            t_ps, d_cr, d_ml = t_ps[mask], d[crater_key][mask], d[melt_key][mask]
        else:
            d_cr, d_ml = d[crater_key], d[melt_key]
        ln1, = ax.plot(t_ps, d_cr, color=st["color"], ls="-", lw=1.8,
                         marker=st["marker"], ms=3, markevery=max(1, len(t_ps) // 12))
        ln2, = ax.plot(t_ps, d_ml, color=st["color"], ls="--", lw=1.0, alpha=0.55)
        legend_crater.append(ln1)
        legend_melt.append(ln2)

    if ep5_markers_ps:
        t1, t2 = ep5_markers_ps
        ax.axvline(t1, color="#d62728", ls=":", lw=1.3, alpha=0.85)
        ax.axvline(t2, color="#9467bd", ls=":", lw=1.3, alpha=0.85)
        ax.text(t1, 0.97, f"Tmax peak\n{t1:.2f} ps", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7, color="#d62728")
        ax.text(t2, 0.82, f"max crater d\n{t2:.2f} ps", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7, color="#9467bd")

    ax.set_xlabel("time (ps)")
    ax.set_ylabel(ylabel)
    title = f"{V26_BANNER}\n{ylabel} vs time"
    if title_suffix:
        title += f"\n{title_suffix}"
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(t_xlim_ps if t_xlim_ps else (0, None))
    ax.set_ylim(bottom=0)
    labels = [CASE_STYLE[c]["label"] for c in (EARLY_PLOT_CASES if t_xlim_ps else DEFAULT_CASE_ORDER)
              if c in case_metrics]
    ax.legend(legend_crater, [f"{lb} crater/vapor" for lb in labels],
              fontsize=6, loc="upper left", framealpha=0.92)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[OK] {out_path}")


def _ep5_reference_times_ps(metrics: dict[str, np.ndarray]) -> tuple[float, float]:
    k_tmax = int(np.argmax(metrics["Tmax_K"]))
    k_maxd = int(np.argmax(metrics["crater_depth_um"]))
    return float(metrics["time_ps"][k_tmax]), float(metrics["time_ps"][k_maxd])


def plot_tfield_panel(
    ax,
    tprint_path: Path,
    mesh: V17MeshGrid,
    metrics: dict[str, np.ndarray],
    frame_mode: str,
    T_melt: float,
    T_vap: float,
    norm: PowerNorm,
) -> object:
    """Draw one mirrored T-field panel; return mappable."""
    times_ns, grids = load_tprint_grids(tprint_path, mesh)
    km = _frame_index_from_metrics(metrics, frame_mode)
    target_ns = float(metrics["time_ns"][km])
    k = int(np.argmin(np.abs(times_ns - target_ns)))
    Tg = grids[k]
    crater_d = float(metrics["crater_depth_um"][km])
    crater_r = float(metrics["crater_radius_um"][km])
    t_ps = float(metrics["time_ps"][km])
    Tmax_K = float(metrics["Tmax_K"][km])

    r_full, T_full = mirror_field(mesh.r_um, Tg)
    z_grid, r_grid = np.meshgrid(mesh.z_um, r_full, indexing="ij")
    im = ax.pcolormesh(r_grid, z_grid, T_full, cmap="inferno", norm=norm, shading="gouraud")
    add_threshold_contours(ax, r_grid, z_grid, T_full, T_melt, T_vap, float(Tg.max()))
    ax.set_xlim(*ZOOM_R_UM)
    ax.set_ylim(*ZOOM_Z_SURFACE_UM)
    ax.set_xlabel("r (um)")
    ax.set_ylabel("z (um)")
    label = "Tmax peak frame" if frame_mode == "tmax_peak" else "max crater-depth frame"
    ax.set_title(
        f"{label}\nt={t_ps:.2f} ps, Tmax={Tmax_K:.0f} K\n"
        f"crater d={crater_d:.3f} um, r={crater_r:.2f} um",
        fontsize=9,
    )
    return im


def plot_v3b_driver_summary(
    case_metrics: dict[str, dict[str, np.ndarray]],
    tprint_path: Path,
    mesh: V17MeshGrid,
    T_melt: float,
    T_vap: float,
    out_path: Path,
) -> None:
    """2x2 V3B driver handoff summary for ep_5p0uj."""
    m = case_metrics[MAIN_CASE]
    t_tmax, t_maxd = _ep5_reference_times_ps(m)

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 10.0))

    for ax, crater_key, melt_key, ylab in (
        (axes[0, 0], "crater_depth_um", "melt_depth_um", "depth (um)"),
        (axes[0, 1], "crater_radius_um", "melt_radius_um", "radius (um)"),
    ):
        t_ps = m["time_ns"] * 1e3
        mask = t_ps <= EARLY_TIME_PS_MAX
        st = CASE_STYLE[MAIN_CASE]
        ax.plot(t_ps[mask], m[crater_key][mask], color=st["color"], ls="-", lw=2.0, label="crater/vapor")
        ax.plot(t_ps[mask], m[melt_key][mask], color=st["color"], ls="--", lw=1.2, alpha=0.7, label="melt")
        ax.axvline(t_tmax, color="#d62728", ls=":", lw=1.4)
        ax.axvline(t_maxd, color="#9467bd", ls=":", lw=1.4)
        ax.set_xlim(0, EARLY_TIME_PS_MAX)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("time (ps)")
        ax.set_ylabel(ylab)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc="upper right")
        ax.set_title(f"{ylab} -- early-time zoom (0--120 ps)", fontsize=10)

    Tmax_all = float(m["Tmax_K"].max())
    norm = PowerNorm(gamma=0.40, vmin=300.0, vmax=Tmax_all)
    im = plot_tfield_panel(
        axes[1, 0], tprint_path, mesh, m, "tmax_peak", T_melt, T_vap, norm,
    )
    plot_tfield_panel(
        axes[1, 1], tprint_path, mesh, m, "max_crater_depth", T_melt, T_vap, norm,
    )

    fig.subplots_adjust(left=0.07, right=0.88, top=0.93, bottom=0.07, hspace=0.38, wspace=0.28)
    cax = fig.add_axes([0.905, 0.12, 0.018, 0.76])
    fig.colorbar(im, cax=cax, label="T (K)")

    fig.suptitle(
        f"{V26_BANNER}\n{MAIN_CASE} V3B driver summary  "
        f"(Tmax peak {t_tmax:.2f} ps | max crater d {t_maxd:.2f} ps)",
        fontsize=11, y=0.98,
    )
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_final_threshold_maps(
    cases: list[str],
    tprint_paths: dict[str, Path],
    case_metrics: dict[str, dict[str, np.ndarray]],
    mesh: V17MeshGrid,
    T_melt: float,
    T_vap: float,
    out_path: Path,
) -> None:
    n = len(cases)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.8 * nrows), squeeze=False)
    Tmax_global = 300.0
    snapshots = {}
    for case in cases:
        m = case_metrics[case]
        k_peak = int(np.argmax(m["Tmax_K"]))
        t_peak = float(m["time_ns"][k_peak])
        times_ns, grids = load_tprint_grids(tprint_paths[case], mesh)
        k = int(np.argmin(np.abs(times_ns - t_peak)))
        Tg = grids[k]
        Tmax_global = max(Tmax_global, float(Tg.max()))
        snapshots[case] = (float(times_ns[k]), Tg)

    norm = PowerNorm(gamma=0.45, vmin=300.0, vmax=Tmax_global)
    im = None
    for ax, case in zip(axes.ravel(), cases):
        t_ns, Tg = snapshots[case]
        r_full, T_full = mirror_field(mesh.r_um, Tg)
        z_grid, r_grid = np.meshgrid(mesh.z_um, r_full, indexing="ij")
        im = ax.pcolormesh(r_grid, z_grid, T_full, cmap="inferno", norm=norm, shading="gouraud")
        handles = add_threshold_contours(ax, r_grid, z_grid, T_full, T_melt, T_vap, float(Tg.max()))
        ax.set_xlim(*ZOOM_R_UM)
        ax.set_ylim(*ZOOM_Z_UM)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        m = case_metrics[case]
        ax.set_title(
            f"{CASE_STYLE[case]['label']}\n"
            f"peak t={t_ns * 1e3:.2f} ps, Tmax={float(Tg.max()):.0f} K",
            fontsize=9,
        )
        ax.text(
            0.02, 0.03,
            f"max melt: d={m['melt_depth_um'].max():.2f} um, r={m['melt_radius_um'].max():.2f} um\n"
            f"max vapor: d={m['vapor_depth_um'].max():.2f} um, r={m['vapor_radius_um'].max():.2f} um",
            transform=ax.transAxes, fontsize=7, color="white", va="bottom",
            bbox=dict(facecolor="black", alpha=0.55, boxstyle="round,pad=0.25"),
        )
        if handles:
            ax.legend(handles=handles, loc="upper right", fontsize=6, framealpha=0.85)
    for ax in axes.ravel()[len(cases):]:
        ax.axis("off")
    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, label="T (K)")
    fig.suptitle(
        f"{V26_BANNER}\nMirrored threshold maps at peak time (top-centre view)",
        y=1.02, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def add_colorbar_right(fig, im, axes, *, label: str = "T (K)") -> None:
    """Dedicated colorbar axis on the far right -- never overlaps subplot grid."""
    fig.subplots_adjust(left=0.06, right=0.88, top=0.90, bottom=0.06, hspace=0.42, wspace=0.28)
    cax = fig.add_axes([0.905, 0.12, 0.018, 0.76])
    cbar = fig.colorbar(im, cax=cax, label=label)
    cbar.ax.tick_params(labelsize=8)


def render_evolution_panels(
    axes_flat,
    times_ns: np.ndarray,
    grids: list[np.ndarray],
    metrics: dict[str, np.ndarray],
    mesh: V17MeshGrid,
    T_melt: float,
    T_vap: float,
    *,
    z_ylim: tuple[float, float],
    target_times_ps: list[float],
) -> object:
    """Fill a flat axes array with mirrored T snapshots; return last mappable."""
    Tmax_all = max(float(g.max()) for g in grids)
    norm = PowerNorm(gamma=0.45, vmin=300.0, vmax=Tmax_all)
    im = None
    for ax, t_ps in zip(axes_flat, target_times_ps):
        target_ns = t_ps * 1e-3
        k, actual_ns, _ = pick_nearest_time_index(times_ns, target_ns, case_name=MAIN_CASE)
        Tg = grids[k]
        r_full, T_full = mirror_field(mesh.r_um, Tg)
        z_grid, r_grid = np.meshgrid(mesh.z_um, r_full, indexing="ij")
        im = ax.pcolormesh(r_grid, z_grid, T_full, cmap="inferno", norm=norm, shading="gouraud")
        add_threshold_contours(ax, r_grid, z_grid, T_full, T_melt, T_vap, float(Tg.max()))
        ax.set_xlim(*ZOOM_R_UM)
        ax.set_ylim(*z_ylim)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        km = int(np.argmin(np.abs(metrics["time_ns"] - actual_ns)))
        ax.set_title(
            f"target={t_ps:g} ps, actual={actual_ns * 1e3:.2f} ps\n"
            f"Tmax={metrics['Tmax_K'][km]:.0f} K, "
            f"crater d={metrics['crater_depth_um'][km]:.3f} um, "
            f"r={metrics['crater_radius_um'][km]:.2f} um",
            fontsize=8,
        )
    return im


def plot_ep5_evolution(
    tprint_path: Path,
    metrics: dict[str, np.ndarray],
    mesh: V17MeshGrid,
    T_melt: float,
    T_vap: float,
    out_path: Path,
) -> None:
    times_ns, grids = load_tprint_grids(tprint_path, mesh)
    n = len(EVOLUTION_TARGETS_PS)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.5, 4.2 * nrows), squeeze=False)
    im = render_evolution_panels(
        axes.ravel(), times_ns, grids, metrics, mesh, T_melt, T_vap,
        z_ylim=ZOOM_Z_UM, target_times_ps=EVOLUTION_TARGETS_PS,
    )
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    if im is not None:
        add_colorbar_right(fig, im, axes)
    fig.suptitle(
        f"{V26_BANNER}\n{MAIN_CASE} crater evolution (mirrored, top-centre)",
        y=0.98, fontsize=11,
    )
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_ep5uj_crater_evolution_zoomz(
    tprint_path: Path,
    metrics: dict[str, np.ndarray],
    mesh: V17MeshGrid,
    T_melt: float,
    T_vap: float,
    out_path: Path,
) -> None:
    """Surface-near vertical zoom so shallow craters (~40 nm) are visible."""
    times_ns, grids = load_tprint_grids(tprint_path, mesh)
    n = len(EVOLUTION_TARGETS_PS)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.5, 4.2 * nrows), squeeze=False)
    im = render_evolution_panels(
        axes.ravel(), times_ns, grids, metrics, mesh, T_melt, T_vap,
        z_ylim=ZOOM_Z_SURFACE_UM, target_times_ps=EVOLUTION_TARGETS_PS,
    )
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    if im is not None:
        add_colorbar_right(fig, im, axes)
    z0, z1 = ZOOM_Z_SURFACE_UM
    fig.suptitle(
        f"{V26_BANNER}\n{MAIN_CASE} crater evolution -- surface zoom "
        f"(z = {z0:.3f}--{z1:.3f} um, dz window = {(z1 - z0) * 1000:.0f} nm)",
        y=0.98, fontsize=11,
    )
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"[OK] {out_path}")


def _frame_index_from_metrics(metrics: dict[str, np.ndarray], mode: str) -> int:
    """Return metrics row index for 'tmax_peak' or 'max_crater_depth' frame."""
    if mode == "tmax_peak":
        return int(np.argmax(metrics["Tmax_K"]))
    if mode == "max_crater_depth":
        return int(np.argmax(metrics["crater_depth_um"]))
    raise ValueError(f"unknown frame mode: {mode}")


def _max_crater_over_time_note(metrics: dict[str, np.ndarray]) -> str:
    """One-line summary of global max crater depth (for Tmax-peak figure footnote)."""
    km = int(np.argmax(metrics["crater_depth_um"]))
    t_ps = float(metrics["time_ps"][km])
    d_um = float(metrics["crater_depth_um"][km])
    return f"max crater depth over time: d={d_um:.3f} um at ~{t_ps:.2f} ps"


def plot_crater_profile_zoom_frame(
    tprint_path: Path,
    metrics: dict[str, np.ndarray],
    mesh: V17MeshGrid,
    T_melt: float,
    T_vap: float,
    out_path: Path,
    *,
    frame_mode: str,
    frame_label: str,
    extra_footnote: str | None = None,
) -> None:
    """
    Zoomed T field + equivalent crater profile for one metrics CSV frame.

    frame_mode: 'tmax_peak' | 'max_crater_depth'
    """
    times_ns, grids = load_tprint_grids(tprint_path, mesh)
    km = _frame_index_from_metrics(metrics, frame_mode)
    target_ns = float(metrics["time_ns"][km])
    k = int(np.argmin(np.abs(times_ns - target_ns)))
    actual_ns = float(times_ns[k])
    Tg = grids[k]

    crater_d = float(metrics["crater_depth_um"][km])
    crater_r = float(metrics["crater_radius_um"][km])
    Tmax_K = float(metrics["Tmax_K"][km])
    t_frame_ps = float(metrics["time_ps"][km])

    fig, (ax_field, ax_prof) = plt.subplots(
        1, 2, figsize=(13.0, 5.2),
        gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.32},
    )

    r_full, T_full = mirror_field(mesh.r_um, Tg)
    z_grid, r_grid = np.meshgrid(mesh.z_um, r_full, indexing="ij")
    norm = PowerNorm(gamma=0.40, vmin=300.0, vmax=float(Tg.max()))
    im = ax_field.pcolormesh(
        r_grid, z_grid, T_full, cmap="inferno", norm=norm, shading="gouraud",
    )
    add_threshold_contours(ax_field, r_grid, z_grid, T_full, T_melt, T_vap, float(Tg.max()))
    ax_field.set_xlim(*ZOOM_R_UM)
    ax_field.set_ylim(*ZOOM_Z_SURFACE_UM)
    ax_field.set_xlabel("r (um)")
    ax_field.set_ylabel("z (um)")
    ax_field.set_title(
        f"T field -- {frame_label}\n"
        f"actual t={t_frame_ps:.2f} ps, Tmax={Tmax_K:.0f} K",
        fontsize=10,
    )
    ax_field.text(
        0.02, 0.04,
        f"crater d={crater_d:.3f} um\n"
        f"crater r={crater_r:.2f} um",
        transform=ax_field.transAxes, fontsize=8, color="white", va="bottom",
        bbox=dict(facecolor="black", alpha=0.6, boxstyle="round,pad=0.25"),
    )
    cax = fig.add_axes([0.46, 0.15, 0.012, 0.72])
    fig.colorbar(im, cax=cax, label="T (K)")

    r_plot = np.linspace(-max(crater_r * 1.4, 25.0), max(crater_r * 1.4, 25.0), 400)
    depth_surface = np.zeros_like(r_plot)
    if crater_d > 0 and crater_r > 0:
        mask = np.abs(r_plot) <= crater_r
        depth_crater = np.zeros_like(r_plot)
        depth_crater[mask] = crater_d * (1.0 - (r_plot[mask] / crater_r) ** 2)
        ax_prof.fill_between(
            r_plot, depth_surface, depth_crater,
            where=mask, color="#d62728", alpha=0.55,
            label=f"equiv. crater (d={crater_d:.3f} um, r={crater_r:.2f} um)",
        )
        ax_prof.plot(r_plot[mask], depth_crater[mask], color="#8b0000", lw=2.0)

    ax_prof.axhline(0.0, color="#333333", lw=1.8, label=f"top surface (z={Z_TOP_UM:.1f} um)")
    if crater_d > 0:
        ax_prof.axhline(crater_d, color="#ff00aa", ls="--", lw=1.2,
                        label=f"frame depth d={crater_d:.3f} um")
    if crater_r > 0:
        ax_prof.axvline(crater_r, color="#888888", ls=":", lw=1.0)
        ax_prof.axvline(-crater_r, color="#888888", ls=":", lw=1.0)
        ax_prof.annotate(
            f"r={crater_r:.2f} um", xy=(crater_r, crater_d * 0.5),
            xytext=(crater_r + 2, min(crater_d * 0.65, DEPTH_PROFILE_MAX_UM * 0.85)),
            fontsize=8,
            arrowprops=dict(arrowstyle="->", color="#555555", lw=0.8),
        )

    j_top = mesh.NZ
    T_top = Tg[j_top, :]
    ax_twin = ax_prof.twinx()
    ax_twin.plot(mesh.r_um, T_top, color="#1f77b4", lw=1.2, alpha=0.85, label="T(r) at top surface")
    ax_twin.axhline(T_vap, color="#ff00aa", ls="--", lw=0.9, alpha=0.7)
    ax_twin.set_ylabel("T at top surface (K)", color="#1f77b4", fontsize=9)
    ax_twin.tick_params(axis="y", labelcolor="#1f77b4", labelsize=8)
    ax_twin.set_xlim(ax_prof.get_xlim())

    ax_prof.set_xlim(-30, 30)
    ax_prof.set_ylim(-0.005, DEPTH_PROFILE_MAX_UM)
    ax_prof.set_xlabel("r (um)")
    ax_prof.set_ylabel("depth below top surface (um)")
    ax_prof.set_title(
        f"Equivalent crater profile -- {frame_label}\n"
        "vertical scale expanded for shallow crater visibility",
        fontsize=10,
    )
    ax_prof.legend(loc="upper left", fontsize=7, framealpha=0.92)
    ax_prof.grid(True, alpha=0.25)

    suptitle = (
        f"{V26_BANNER}\n{MAIN_CASE} {frame_label}  "
        f"(t={t_frame_ps:.2f} ps, Tmax={Tmax_K:.0f} K, "
        f"d={crater_d:.3f} um, r={crater_r:.2f} um)"
    )
    fig.suptitle(suptitle, fontsize=10, y=1.03)

    if extra_footnote:
        fig.text(
            0.5, 0.01, extra_footnote,
            ha="center", va="bottom", fontsize=8.5, color="#333",
            bbox=dict(facecolor="#f5f5f5", edgecolor="#cccccc", boxstyle="round,pad=0.35"),
        )

    fig.subplots_adjust(left=0.07, right=0.88, top=0.86, bottom=0.14 if extra_footnote else 0.12,
                      wspace=0.38)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_peak_crater_profile_zoom(
    tprint_path: Path,
    metrics: dict[str, np.ndarray],
    mesh: V17MeshGrid,
    _summary_row: dict[str, str],
    T_melt: float,
    T_vap: float,
    out_path: Path,
) -> None:
    """Tmax peak frame (not necessarily max crater depth)."""
    km_tmax = _frame_index_from_metrics(metrics, "tmax_peak")
    t_ps = float(metrics["time_ps"][km_tmax])
    d_ps = float(metrics["crater_depth_um"][km_tmax])
    footnote = (
        f"Tmax peak frame: t={t_ps:.2f} ps, crater d={d_ps:.3f} um  |  "
        f"{_max_crater_over_time_note(metrics)}"
    )
    plot_crater_profile_zoom_frame(
        tprint_path, metrics, mesh, T_melt, T_vap, out_path,
        frame_mode="tmax_peak",
        frame_label="Tmax peak frame",
        extra_footnote=footnote,
    )


def plot_max_crater_profile_zoom(
    tprint_path: Path,
    metrics: dict[str, np.ndarray],
    mesh: V17MeshGrid,
    T_melt: float,
    T_vap: float,
    out_path: Path,
) -> None:
    """Frame where crater_depth_um is global maximum over time."""
    plot_crater_profile_zoom_frame(
        tprint_path, metrics, mesh, T_melt, T_vap, out_path,
        frame_mode="max_crater_depth",
        frame_label="max crater-depth frame",
        extra_footnote=None,
    )


def plot_revolved_crater_schematic(
    summary_row: dict[str, str],
    out_path: Path,
) -> None:
    """Post-process schematic: top-view disk + side cross-section + simple 3D surface."""
    depth = float(summary_row["max_crater_depth_um"])
    radius = float(summary_row["max_crater_radius_um"])
    t_peak_ps = float(summary_row["t_peak_ps"])
    T_peak = float(summary_row["Tmax_peak_K"])

    fig = plt.figure(figsize=(13.0, 5.0))
    wafer_r = 100.0

    ax1 = fig.add_subplot(1, 3, 1)
    theta = np.linspace(0, 2 * np.pi, 200)
    ax1.fill(wafer_r * np.cos(theta), wafer_r * np.sin(theta),
             color="#cfd8dc", ec="#546e7a", lw=1.2, label="Si wafer (100 um radius)")
    if radius > 0:
        ax1.fill(radius * np.cos(theta), radius * np.sin(theta),
                 color="#d62728", alpha=0.65, ec="#8b0000", lw=1.5,
                 label=f"equiv. crater r={radius:.2f} um")
    ax1.plot(0, 0, "k+", ms=8)
    ax1.set_aspect("equal")
    ax1.set_xlabel("x (um)")
    ax1.set_ylabel("y (um)")
    ax1.set_title("Top view (revolved)\ncrater at wafer centre")
    ax1.legend(fontsize=7, loc="upper right")
    ax1.grid(True, alpha=0.25)

    ax2 = fig.add_subplot(1, 3, 2)
    z_top = 5.0
    ax2.axhspan(0, z_top, color="#cfd8dc", alpha=0.35)
    ax2.plot([-wafer_r, wafer_r], [z_top, z_top], "k-", lw=1.5, label="top surface")
    if depth > 0 and radius > 0:
        r_cr = np.linspace(-radius, radius, 100)
        z_cr = z_top - depth * (1.0 - (r_cr / radius) ** 2)
        z_cr = np.clip(z_cr, z_top - depth, z_top)
        ax2.fill_between(r_cr, z_cr, z_top, color="#d62728", alpha=0.6,
                         label=f"equiv. crater d={depth:.2f} um")
    ax2.axvline(0, color="gray", ls="--", lw=0.8)
    ax2.set_xlim(-60, 60)
    ax2.set_ylim(4.85, 5.02)
    ax2.set_xlabel("r (um)")
    ax2.set_ylabel("z (um)")
    ax2.set_title("Side cross-section\n(parabolic schematic)")
    ax2.legend(fontsize=7)
    ax2.grid(True, alpha=0.25)

    ax3 = fig.add_subplot(1, 3, 3, projection="3d")
    rr = np.linspace(-wafer_r, wafer_r, 80)
    tt = np.linspace(0, 2 * np.pi, 80)
    RR, TT = np.meshgrid(rr, tt)
    XX, YY = RR * np.cos(TT), RR * np.sin(TT)
    ZZ = np.full_like(XX, z_top)
    if depth > 0 and radius > 0:
        crater_mask = np.sqrt(XX ** 2 + YY ** 2) <= radius
        r_xy = np.sqrt(XX ** 2 + YY ** 2)
        ZZ[crater_mask] = z_top - depth * (1.0 - (r_xy[crater_mask] / radius) ** 2)
    ax3.plot_surface(XX, YY, ZZ, cmap="coolwarm", alpha=0.85, linewidth=0, antialiased=True)
    ax3.set_xlabel("x (um)")
    ax3.set_ylabel("y (um)")
    ax3.set_zlabel("z (um)")
    ax3.set_title("3D revolved schematic\n(not LS-DYNA 3D model)")
    ax3.view_init(elev=35, azim=-60)

    fig.suptitle(
        f"{V26_BANNER}\n{MAIN_CASE} peak crater schematic  "
        f"(t_peak={t_peak_ps:.1f} ps, Tmax={T_peak:.0f} K, "
        f"d={depth:.2f} um, r={radius:.2f} um)",
        fontsize=10, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--toml", type=Path,
                    default=project_root / "config" / "v17_30ps_local_mesh.toml")
    ap.add_argument("--registry", type=Path,
                    default=project_root / "models" / "v17_30ps_local" / "v17_case_registry.csv")
    ap.add_argument("--mesh-nodes", type=Path,
                    default=project_root / "models" / "v17_30ps_local" / "v17_mesh_nodes.csv")
    ap.add_argument("--v26-root", type=Path,
                    default=project_root / "results" / "v26_30ps_threshold_ablation")
    ap.add_argument("--v17-root", type=Path,
                    default=project_root / "results" / "v17_30ps_local")
    args = ap.parse_args()

    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)
    T_melt = float(cfg["material"]["T_melt_K"])
    T_vap = float(cfg["material"]["T_vap_K"])

    summary_csv = args.v26_root / "v26_case_summary.csv"
    if not summary_csv.is_file():
        print(f"[ERROR] run extract first: {summary_csv}", file=sys.stderr)
        return 2

    with args.registry.open("r", encoding="utf-8") as fh:
        registry = {row["name"]: row for row in csv.DictReader(fh)}
    mesh = load_mesh_grid(args.mesh_nodes, registry[MAIN_CASE])

    summary_rows = read_summary_csv(summary_csv)
    case_metrics: dict[str, dict[str, np.ndarray]] = {}
    tprint_paths: dict[str, Path] = {}
    available: list[str] = []

    for case in DEFAULT_CASE_ORDER:
        mcsv = args.v26_root / case / "v26_threshold_metrics.csv"
        tp = args.v17_root / case / "tprint"
        if not mcsv.is_file():
            print(f"  [SKIP] {case}: missing metrics CSV", file=sys.stderr)
            continue
        if not tp.is_file():
            print(f"  [SKIP] {case}: missing tprint", file=sys.stderr)
            continue
        case_metrics[case] = read_metrics_csv(mcsv)
        tprint_paths[case] = tp
        available.append(case)

    if not available:
        print("[ERROR] no usable cases for plotting", file=sys.stderr)
        return 3

    fig_dir = args.v26_root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_case_comparison_table(summary_rows, fig_dir / "v26_case_comparison_table.png")
    plot_metric_vs_time(
        case_metrics, "crater_depth_um", "melt_depth_um", "depth (um)",
        fig_dir / "v26_crater_depth_vs_time.png",
    )
    plot_metric_vs_time(
        case_metrics, "crater_radius_um", "melt_radius_um", "radius (um)",
        fig_dir / "v26_crater_radius_vs_time.png",
    )
    if MAIN_CASE in case_metrics:
        t_tmax, t_maxd = _ep5_reference_times_ps(case_metrics[MAIN_CASE])
        plot_metric_vs_time(
            case_metrics, "crater_depth_um", "melt_depth_um", "depth (um)",
            fig_dir / "v26_crater_depth_vs_time_early.png",
            t_xlim_ps=(0.0, EARLY_TIME_PS_MAX),
            title_suffix="early-time zoom (0--120 ps)",
            ep5_markers_ps=(t_tmax, t_maxd),
        )
        plot_metric_vs_time(
            case_metrics, "crater_radius_um", "melt_radius_um", "radius (um)",
            fig_dir / "v26_crater_radius_vs_time_early.png",
            t_xlim_ps=(0.0, EARLY_TIME_PS_MAX),
            title_suffix="early-time zoom (0--120 ps)",
            ep5_markers_ps=(t_tmax, t_maxd),
        )
    plot_final_threshold_maps(
        [c for c in DEFAULT_CASE_ORDER if c in available],
        tprint_paths, case_metrics, mesh, T_melt, T_vap,
        fig_dir / "v26_final_threshold_maps.png",
    )

    if MAIN_CASE in available:
        ep5_summary = next(
            r for r in summary_rows if r["case"] == MAIN_CASE and r.get("status") == "OK"
        )
        plot_ep5_evolution(
            tprint_paths[MAIN_CASE], case_metrics[MAIN_CASE], mesh,
            T_melt, T_vap, fig_dir / "v26_ep5uj_crater_evolution.png",
        )
        plot_ep5uj_crater_evolution_zoomz(
            tprint_paths[MAIN_CASE], case_metrics[MAIN_CASE], mesh,
            T_melt, T_vap, fig_dir / "v26_ep5uj_crater_evolution_zoomz.png",
        )
        plot_peak_crater_profile_zoom(
            tprint_paths[MAIN_CASE], case_metrics[MAIN_CASE], mesh,
            ep5_summary, T_melt, T_vap,
            fig_dir / "v26_ep5uj_peak_crater_profile_zoom.png",
        )
        plot_max_crater_profile_zoom(
            tprint_paths[MAIN_CASE], case_metrics[MAIN_CASE], mesh,
            T_melt, T_vap,
            fig_dir / "v26_ep5uj_max_crater_profile_zoom.png",
        )
        plot_v3b_driver_summary(
            case_metrics, tprint_paths[MAIN_CASE], mesh,
            T_melt, T_vap, fig_dir / "v26_ep5p0uj_v3b_driver_summary.png",
        )
        plot_revolved_crater_schematic(
            ep5_summary, fig_dir / "v26_ep5uj_revolved_crater_schematic.png",
        )
    else:
        print(f"  [WARN] {MAIN_CASE} not available; skipping driver figures", file=sys.stderr)

    print(f"\n[DONE] figures in {fig_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
