"""
plot_v5b_unified_ablation_sequence.py
=====================================

Unified early ablation + plume diffusion sequence on a single fixed coordinate
frame (r, z).  One continuous simulation scene -- not separate plume/substrate
panels stitched together.

Run:
    python scripts\\build_v3b_plume_model.py
    python scripts\\plot_v5b_unified_ablation_sequence.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

_THIS = Path(__file__).resolve().parent
if str(_THIS) not in sys.path:
    sys.path.insert(0, str(_THIS))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
from matplotlib import font_manager  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

from plot_v3b_plume_shock import (  # noqa: E402
    generate_ejecta_particles,
    metric_R0,
    read_metrics,
    read_t0_ps,
)

# ---------------------------------------------------------------------------
# Fixed layout constants (never change per frame)
# ---------------------------------------------------------------------------
R_LIM = (-60.0, 60.0)
Z_LIM = (-8.0, 80.0)
MAIN_TIMES_PS = [0, 10, 30, 50, 100, 200, 500]
LATE_TIMES_PS = [1000, 2000, 5000]
LATE_Z_LIM = (-8.0, 200.0)

CRATER_DEPTH_EXAG = 60.0  # fixed visual exag, same all frames
T_SUB_VMIN = 300.0
T_SUB_VMAX = 5600.0
PLUME_CMAP = "magma"
SUB_CMAP = "copper"
BG_COLOR = "#f4f4f4"
MAIN_TITLE = "30 ps 激光烧蚀下硅片早期羽流扩散演化（ep_5p0uj, early stage）"
LATE_TITLE = "30 ps 激光烧蚀硅片羽流扩散（ep_5p0uj, late supplement）"


def configure_fonts() -> None:
    for name in ("Microsoft YaHei", "SimHei", "SimSun", "Noto Sans CJK SC", "DejaVu Sans"):
        if name in {f.name for f in font_manager.fontManager.ttflist}:
            plt.rcParams["font.sans-serif"] = [name]
            break
    plt.rcParams["axes.unicode_minus"] = False


configure_fonts()


@dataclass(frozen=True)
class FrameGeom:
    t_ps: float
    R_cr_um: float
    d_cr_um: float
    d_cr_draw_um: float


def load_driver(root: Path) -> dict:
    for rel in (
        "results/v26_30ps_threshold_ablation/v3b_driver/v3b_driver_ep5p0uj.json",
        "results/v3b_30ps_plume/v3b_driver_used.json",
    ):
        p = root / rel
        if p.is_file():
            return json.loads(p.read_text(encoding="utf-8"))
    return {}


def load_v26_crater(root: Path) -> dict[str, np.ndarray] | None:
    p = root / "results/v26_30ps_threshold_ablation/ep_5p0uj/v26_threshold_metrics.csv"
    if not p.is_file():
        return None
    with p.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return {
        "time_ps": np.array([float(r["time_ps"]) for r in rows]),
        "crater_radius_um": np.array([float(r["crater_radius_um"]) for r in rows]),
        "crater_depth_um": np.array([float(r["crater_depth_um"]) for r in rows]),
        "Tmax_K": np.array([float(r["Tmax_K"]) for r in rows]),
    }


def nearest_index(times: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(times - target)))


def interp_crater(v26: dict[str, np.ndarray] | None, t_ps: float) -> tuple[float, float, float]:
    if v26 is None:
        return 0.0, 0.0, 300.0
    t = v26["time_ps"]
    R = float(np.interp(t_ps, t, v26["crater_radius_um"], left=0.0, right=v26["crater_radius_um"][-1]))
    d = float(np.interp(t_ps, t, v26["crater_depth_um"], left=0.0, right=v26["crater_depth_um"][-1]))
    Tmax = float(np.interp(t_ps, t, v26["Tmax_K"], left=300.0, right=v26["Tmax_K"][-1]))
    return R, d, Tmax


def frame_geometry(v26: dict | None, t_ps: float) -> FrameGeom:
    R_cr, d_cr, _ = interp_crater(v26, t_ps)
    d_draw = d_cr * CRATER_DEPTH_EXAG if d_cr > 1e-9 else 0.0
    if d_cr > 1e-9:
        d_draw = max(d_draw, 0.12)  # fixed floor so shallow crater is visible in [-8,80] frame
    return FrameGeom(t_ps=t_ps, R_cr_um=R_cr, d_cr_um=d_cr, d_cr_draw_um=d_draw)


def surface_z_profile(r: np.ndarray, geom: FrameGeom) -> np.ndarray:
    """Original surface z=0 with axisymmetric parabolic crater bowl (z <= 0)."""
    if geom.R_cr_um <= 1e-9 or geom.d_cr_draw_um <= 0:
        return np.zeros_like(r)
    rr = np.abs(r)
    z_s = np.zeros_like(r)
    inside = rr <= geom.R_cr_um
    z_s[inside] = -geom.d_cr_draw_um * (1.0 - (rr[inside] / geom.R_cr_um) ** 2)
    return z_s


def substrate_temperature(
    R: np.ndarray,
    Z: np.ndarray,
    geom: FrameGeom,
    v26: dict | None,
) -> np.ndarray:
    """Temperature proxy inside silicon (below local surface)."""
    r = R
    z_surf = surface_z_profile(r, geom)
    inside = Z <= z_surf
    if not np.any(inside):
        return np.full_like(R, np.nan)

    _, _, T_peak = interp_crater(v26, geom.t_ps)
    depth = np.clip(z_surf - Z, 0.0, None)
    rr = np.abs(r)
    R_m = max(float(interp_crater(v26, min(geom.t_ps, 30.0))[0]), 8.0) if v26 else 15.0
    T = T_SUB_VMIN + (T_peak - T_SUB_VMIN) * np.exp(-rr ** 2 / max(R_m, 1.0) ** 2)
    T *= np.exp(-depth / max(0.35 + 0.004 * geom.t_ps, 0.15))
    out = np.full_like(R, np.nan)
    out[inside] = T[inside]
    return out


def plume_density_coupled(
    R: np.ndarray,
    Z: np.ndarray,
    geom: FrameGeom,
    M_row: dict[str, float],
    t0_ps: float,
    R0_launch: float,
) -> np.ndarray:
    """
    Plume density anchored to crater mouth / surface (z >= local surface).
    Uses V3B metrics but couples emission to crater radius and launch fraction.
    """
    r = R
    z_surf = surface_z_profile(r, geom)
    above = Z >= z_surf
    if not np.any(above):
        return np.zeros_like(R)

    t_ps = geom.t_ps
    launch_frac = 0.0 if t_ps < t0_ps - 1e-6 else min(1.0, ((t_ps - t0_ps) / max(t0_ps, 1.0)) ** 0.55 + 0.15)
    if t_ps < t0_ps:
        # Pre-launch: weak vapor wisps from heating (scales with crater growth).
        launch_frac = min(0.35, max(geom.d_cr_um / 0.042, 0.0) * (t_ps / max(t0_ps, 1.0)) ** 1.5)

    n0 = float(M_row["n0_rel"]) * launch_frac
    if n0 <= 0 and launch_frac <= 0:
        return np.zeros_like(R)

    sigma_r = max(float(M_row["sigma_r_um"]), 2.0)
    sigma_z = max(float(M_row["sigma_z_um"]), 2.0)
    z_center = max(float(M_row["z_center_um"]), 1.5)
    R_mouth = max(geom.R_cr_um, R0_launch * 0.35, 1.0)

    rr = np.abs(r)
    zz = Z - z_surf  # height above local surface
    zz = np.clip(zz, 0.0, None)

    # Main plume blob rising from mouth.
    blob = n0 * np.exp(-(rr / sigma_r) ** 2 - ((zz - z_center) / sigma_z) ** 2)
    # Crater-root emission (surface attachment).
    root = (
        launch_frac
        * np.exp(-(rr / (0.85 * R_mouth)) ** 2)
        * np.exp(-(zz / np.maximum(0.35 + 0.08 * zz, 0.25)) ** 2)
    )
    density = blob + 0.55 * n0 * root
    # Mouth collar: extra density at lip.
    lip = 0.25 * n0 * np.exp(-((rr - 0.85 * R_mouth) / max(R_mouth * 0.25, 0.5)) ** 2) * np.exp(-(zz / 1.2) ** 2)
    density = density + lip

    out = np.zeros_like(R)
    out[above] = density[above]
    return out


def silicon_solid_mask(R: np.ndarray, Z: np.ndarray, geom: FrameGeom) -> np.ndarray:
    return Z <= surface_z_profile(R, geom)


def precompute_plume_vmax(
    times_ps: list[float],
    M: dict,
    v26: dict | None,
    t0_ps: float,
    R0: float,
    r: np.ndarray,
    z: np.ndarray,
) -> float:
    Rg, Zg = np.meshgrid(r, z, indexing="xy")
    vmax = 1e-6
    for t in times_ps:
        k = nearest_index(M["time_ps"], t)
        geom = frame_geometry(v26, t)
        row = {key: float(M[key][k]) for key in M}
        d = plume_density_coupled(Rg, Zg, geom, row, t0_ps, R0)
        vmax = max(vmax, float(d.max()))
    return vmax


def draw_crater_silicon_fill(ax, geom: FrameGeom, rlim: tuple[float, float], zlim: tuple[float, float]) -> None:
    """Fill silicon solid below crater surface (subtle base tone)."""
    r_line = np.linspace(rlim[0], rlim[1], 400)
    z_s = surface_z_profile(r_line, geom)
    verts = (
        list(zip(r_line, z_s))
        + [(rlim[1], zlim[0]), (rlim[0], zlim[0])]
    )
    ax.add_patch(Polygon(verts, closed=True, facecolor="#8d8d8d", edgecolor="none", alpha=0.35, zorder=1))


def draw_surface_line(ax, geom: FrameGeom, rlim: tuple[float, float]) -> None:
    r_line = np.linspace(rlim[0], rlim[1], 400)
    z_s = surface_z_profile(r_line, geom)
    ax.plot(r_line, z_s, color="#2a2a2a", lw=1.0, zorder=4)


def draw_plume_front(ax, geom: FrameGeom, R_plume: float, t0_ps: float) -> None:
    if geom.t_ps < t0_ps - 1e-6 or R_plume <= 0:
        return
    th = np.linspace(0, np.pi, 100)
    z_surf_r = 0.0
    ax.plot(R_plume * np.cos(th), z_surf_r + R_plume * np.sin(th) * 0.55,
            color="#00acc1", ls="--", lw=0.55, alpha=0.55, zorder=6)


def draw_shock_front(ax, geom: FrameGeom, R_shock: float, R_plume: float, t0_ps: float) -> None:
    if geom.t_ps < t0_ps - 1e-6 or R_shock <= R_plume + 0.5:
        return
    th = np.linspace(0, np.pi, 100)
    ax.plot(R_shock * np.cos(th), R_shock * np.sin(th) * 0.55,
            color="#ef5350", ls="-", lw=0.45, alpha=0.45, zorder=6)


def draw_ejecta_sparse(
    ax,
    geom: FrameGeom,
    M_row: dict,
    t0_ps: float,
    seed: int,
    frame_idx: int,
) -> None:
    if geom.t_ps < t0_ps - 1e-6:
        return
    eh = float(M_row["ejecta_height_um"])
    er = float(M_row["ejecta_radius_um"])
    if eh <= 0.2:
        return
    dt_ns = max((geom.t_ps - t0_ps) * 1e-3, 0.0)
    strength = (1.0 - np.exp(-dt_ns / 0.15)) * np.exp(-dt_ns / 2.5)
    if strength < 0.08:
        return
    n_particles = max(int(60 * strength), 8)
    rng = np.random.default_rng(seed + frame_idx)
    x, z, _ = generate_ejecta_particles(
        n_particles, max(geom.R_cr_um, 1.0), eh * 0.35, er, 0.0, rng, strength * 0.6,
    )
    if x.size:
        ax.scatter(x, z, s=2.2, c="#ffb74d", alpha=0.28, linewidths=0, zorder=5)


def draw_frame(
    ax,
    t_tgt_ps: float,
    M: dict,
    v26: dict | None,
    t0_ps: float,
    R0: float,
    r: np.ndarray,
    z: np.ndarray,
    plume_norm: Normalize,
    sub_norm: Normalize,
    rlim: tuple[float, float],
    zlim: tuple[float, float],
    frame_idx: int,
    ejecta_seed: int,
) -> tuple[object, object]:
    k = nearest_index(M["time_ps"], t_tgt_ps)
    geom = frame_geometry(v26, t_tgt_ps)
    row = {key: float(M[key][k]) for key in M}
    Rg, Zg = np.meshgrid(r, z, indexing="xy")

    ax.set_facecolor(BG_COLOR)

    # 1. Silicon base fill.
    draw_crater_silicon_fill(ax, geom, rlim, zlim)

    # 2. Substrate internal field.
    T_sub = substrate_temperature(Rg, Zg, geom, v26)
    sub_im = ax.pcolormesh(r, z, T_sub, cmap=SUB_CMAP, norm=sub_norm, shading="gouraud", zorder=2)

    # 3. Plume above surface (coupled).
    dens = plume_density_coupled(Rg, Zg, geom, row, t0_ps, R0)
    plume_im = ax.pcolormesh(
        r, z, np.ma.masked_where(dens <= 1e-8, dens),
        cmap=PLUME_CMAP, norm=plume_norm, shading="gouraud", alpha=0.88, zorder=3,
    )

    # 4. Surface / crater line.
    draw_surface_line(ax, geom, rlim)

    # 5. Sparse ejecta.
    draw_ejecta_sparse(ax, geom, row, t0_ps, ejecta_seed, frame_idx)

    # 6. Thin fronts.
    draw_plume_front(ax, geom, float(row["R_plume_um"]), t0_ps)
    draw_shock_front(ax, geom, float(row["R_shock_um"]), float(row["R_plume_um"]), t0_ps)

    ax.set_xlim(rlim)
    ax.set_ylim(zlim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{int(t_tgt_ps)} ps", fontsize=9, pad=3)
    return sub_im, plume_im


def style_axes(ax, show_ylabel: bool = False) -> None:
    ax.tick_params(labelsize=7)
    if show_ylabel:
        ax.set_ylabel("z / um", fontsize=8)
    ax.set_xlabel("signed r / um", fontsize=8)


def build_legend_panel(ax, t0_ps: float) -> None:
    ax.set_facecolor(BG_COLOR)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    lines = [
        Line2D([0], [0], color="#8d8d8d", lw=6, alpha=0.5, label="Si substrate"),
        Line2D([0], [0], color="#2a2a2a", lw=1.2, label="Surface / crater"),
        Line2D([0], [0], color="#d84315", lw=6, alpha=0.8, label="Plume density proxy"),
        Line2D([0], [0], color="#00acc1", ls="--", lw=0.8, label="Plume front"),
        Line2D([0], [0], color="#ef5350", lw=0.8, label="Shock front"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#ffb74d", markersize=4,
               alpha=0.5, linestyle="None", label="Ejecta proxy (sparse)"),
    ]
    ax.legend(handles=lines, loc="center", fontsize=7, framealpha=0.95, title="Layers")
    txt = (
        f"Unified frame: r=[{R_LIM[0]:.0f},{R_LIM[1]:.0f}] um\n"
        f"z=[{Z_LIM[0]:.0f},{Z_LIM[1]:.0f}] um, z=0 surface\n"
        f"Plume launch t0={t0_ps:.1f} ps\n"
        f"Crater from V2.6; plume from V3B\n"
        f"Parametric proxy (not MD/SPH)"
    )
    ax.text(0.5, 0.08, txt, transform=ax.transAxes, ha="center", va="bottom", fontsize=7, color="0.25")


def plot_sequence(
    times_ps: list[float],
    out_path: Path,
    title: str,
    M: dict,
    v26: dict | None,
    t0_ps: float,
    R0: float,
    rlim: tuple[float, float],
    zlim: tuple[float, float],
    layout: str,
    ejecta_seed: int,
) -> None:
    r = np.linspace(rlim[0], rlim[1], 360)
    z = np.linspace(zlim[0], zlim[1], 240)
    plume_vmax = precompute_plume_vmax(times_ps, M, v26, t0_ps, R0, r, z)
    plume_norm = Normalize(vmin=0.0, vmax=plume_vmax)
    sub_norm = Normalize(vmin=T_SUB_VMIN, vmax=T_SUB_VMAX)

    if layout == "main":
        fig, axes = plt.subplots(2, 4, figsize=(14.0, 7.2), dpi=300)
        panel_axes = list(axes.ravel()[:7])
        legend_ax = axes.ravel()[7]
    else:
        n = len(times_ps)
        fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 4.6), dpi=300)
        panel_axes = list(np.atleast_1d(axes).ravel())
        legend_ax = None

    sub_im = plume_im = None
    ncols_row = 4 if layout == "main" else len(times_ps)
    for i, (ax, t_ps) in enumerate(zip(panel_axes, times_ps)):
        sub_im, plume_im = draw_frame(
            ax, t_ps, M, v26, t0_ps, R0, r, z, plume_norm, sub_norm,
            rlim, zlim, i, ejecta_seed,
        )
        style_axes(ax, show_ylabel=(i % ncols_row == 0))

    if legend_ax is not None:
        build_legend_panel(legend_ax, t0_ps)

    if layout == "main":
        fig.subplots_adjust(left=0.05, right=0.88, top=0.90, bottom=0.10, hspace=0.22, wspace=0.12)
        if sub_im is not None:
            cax1 = fig.add_axes([0.895, 0.58, 0.012, 0.28])
            fig.colorbar(sub_im, cax=cax1, label="Substrate T proxy (K)")
        if plume_im is not None:
            cax2 = fig.add_axes([0.895, 0.16, 0.012, 0.28])
            fig.colorbar(plume_im, cax=cax2, label="Plume density proxy")
    else:
        fig.subplots_adjust(left=0.06, right=0.90, top=0.88, bottom=0.14, wspace=0.18)
        if sub_im is not None:
            cax1 = fig.add_axes([0.915, 0.58, 0.012, 0.28])
            fig.colorbar(sub_im, cax=cax1, label="Substrate T proxy (K)")
        if plume_im is not None:
            cax2 = fig.add_axes([0.915, 0.16, 0.012, 0.28])
            fig.colorbar(plume_im, cax=cax2, label="Plume density proxy")

    fig.suptitle(title, fontsize=12, y=0.98)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_summary(out_path: Path, M: dict, v26: dict | None, t0_ps: float) -> None:
    """Small separate summary: crater radius/depth vs time."""
    if v26 is None:
        return
    fig, ax1 = plt.subplots(figsize=(6.5, 3.8), dpi=200)
    ax2 = ax1.twinx()
    t = v26["time_ps"]
    m = t <= 600
    ax1.plot(t[m], v26["crater_radius_um"][m], "b-o", ms=3, lw=1.2, label="crater radius")
    ax2.plot(t[m], v26["crater_depth_um"][m] * CRATER_DEPTH_EXAG, "r-s", ms=3, lw=1.2,
             label=f"crater depth x{CRATER_DEPTH_EXAG:.0f} (draw)")
    for tt in MAIN_TIMES_PS:
        ax1.axvline(tt, color="0.85", lw=0.5)
    ax1.axvline(t0_ps, color="#00acc1", ls="--", lw=0.8, label=f"t0={t0_ps:.1f} ps")
    ax1.set_xlabel("time (ps)")
    ax1.set_ylabel("crater radius (um)")
    ax2.set_ylabel("drawn crater depth (um)")
    ax1.set_title("V2.6 crater evolution (summary, separate from main sequence)")
    ax1.grid(True, alpha=0.3)
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [ln.get_label() for ln in lines], fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()

    csv_path = root / "results/v3b_30ps_plume/v3b_plume_shock_metrics.csv"
    t0_path = root / "results/v3b_30ps_plume/v3b_t0.txt"
    if not csv_path.is_file():
        print(f"[ERROR] missing {csv_path}; run build_v3b_plume_model.py", file=sys.stderr)
        return 2

    M = read_metrics(csv_path)
    t0_ps = read_t0_ps(t0_path)
    R0 = metric_R0(M)
    v26 = load_v26_crater(root)
    out_dir = args.out_dir or (root / "figures/v5b_reference_style/ep_5p0uj_ejecta_8500ms")
    ejecta_seed = 42

    plot_sequence(
        MAIN_TIMES_PS,
        out_dir / "v5b_early_ablation_plume_sequence.png",
        MAIN_TITLE,
        M, v26, t0_ps, R0, R_LIM, Z_LIM, "main", ejecta_seed,
    )
    plot_sequence(
        LATE_TIMES_PS,
        out_dir / "v5b_late_ablation_plume_sequence.png",
        LATE_TITLE,
        M, v26, t0_ps, R0, R_LIM, LATE_Z_LIM, "late", ejecta_seed,
    )
    plot_summary(out_dir / "v5b_crater_evolution_summary.png", M, v26, t0_ps)

    print(f"\n[DONE] unified V5B sequence -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
