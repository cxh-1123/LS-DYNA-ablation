"""
explain_axisymmetric_center_v16.py
==================================

V1.6 step 3 -- a teaching figure that explains why the heated region
shows up in the *upper-left corner* of LS-PrePost's V1/V1.5 view, and
how the mirror / 3-D revolve operations map it back to the *centre of
the silicon disk*.

Outputs (both PNG at 300 dpi, English labels only):

  results/v16_30ps/figures/v16_axisymmetric_coordinate_explanation.png
      original 1 x 3 full-view triptych (unchanged purpose).

  results/v16_30ps/figures/v16_axisymmetric_coordinate_explanation_zoomed.png
      2 x 3 composite: row 1 = full views, row 2 = zoomed views of the
      axis / top-surface hot zone for clearer presentation.

By default we re-use the V2 vapor_confirm_100ns peak snapshot
(read-only) for an authentic temperature field; if that file is
missing we fall back to a synthetic Gaussian heat spot.  V1/V1.5/V2/
V3A files are NOT modified.

Run from project root:
    .\\.venv\\Scripts\\python.exe scripts\\explain_axisymmetric_center_v16.py
"""

from __future__ import annotations

import argparse
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.ticker import FuncFormatter


_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))


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
    "legend.fontsize": 8,
    "savefig.dpi": 300,
})

# Zoom windows (um) -- fixed for reproducible presentation figures
ZOOM_R_HALF_UM = 50.0          # Panel 1 zoom: r = 0 .. 50 um
ZOOM_R_FULL_UM = 50.0          # Panel 2 zoom: r = -50 .. +50 um
ZOOM_Z_MIN_UM = 180.0          # both section zooms: z = 180 .. 200 um
ZOOM_Z_MAX_UM = 200.0
ZOOM_DISK_HALF_UM = 100.0      # Panel 3 zoom: x/y = -100 .. +100 um


@dataclass
class FieldData:
    r_mm: np.ndarray
    z_mm: np.ndarray
    T_grid: np.ndarray
    source_label: str
    using_real: bool


# =============================================================================
# Data loading
# =============================================================================
def load_v2_peak_snapshot(project_root: Path, case: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from lsdyna_tprint import (
        RADIUS_MM, THICKNESS_MM,
        parse_tprint, reshape_to_grid,
    )

    tprint = project_root / "results" / "v15" / case / "tprint"
    if not tprint.is_file():
        raise FileNotFoundError(f"V1.5 tprint not found for {case}: {tprint}")

    times_ms, T_all = parse_tprint(tprint)
    k_peak = int(np.argmax(T_all.max(axis=1)))
    T_grid = reshape_to_grid(T_all[k_peak])
    r_mm = np.linspace(0.0, RADIUS_MM, T_grid.shape[1])
    z_mm = np.linspace(0.0, THICKNESS_MM, T_grid.shape[0])
    return r_mm, z_mm, T_grid


def synth_field(
    radius_um: float, thickness_um: float,
    w0_um: float, peak_T_K: float, base_T_K: float = 300.0,
    NR: int = 100, NZ: int = 80,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r_um = np.linspace(0.0, radius_um, NR + 1)
    z_um = np.linspace(0.0, thickness_um, NZ + 1)
    R, Z = np.meshgrid(r_um, z_um, indexing="xy")
    decay_um = max(1.0, 0.02 * thickness_um)
    spatial = np.exp(-2.0 * R ** 2 / w0_um ** 2)
    depth = np.exp(-(thickness_um - Z) / decay_um)
    T = base_T_K + (peak_T_K - base_T_K) * spatial * depth
    return r_um / 1000.0, z_um / 1000.0, T


def load_field_data(project_root: Path, cfg: dict) -> FieldData:
    case = cfg["axisymmetric_explanation"]["prefer_v2_case"]
    try:
        r_mm, z_mm, T_grid = load_v2_peak_snapshot(project_root, case)
        return FieldData(
            r_mm, z_mm, T_grid,
            source_label=f"data: V1.5 / V2 case `{case}` peak snapshot",
            using_real=True,
        )
    except (FileNotFoundError, ImportError) as e:
        w0_um = float(cfg["axisymmetric_explanation"]["fallback_synthetic_w0_um"])
        peak_T = float(cfg["axisymmetric_explanation"]["fallback_peak_T_K"])
        r_mm, z_mm, T_grid = synth_field(
            radius_um=500.0, thickness_um=200.0,
            w0_um=w0_um, peak_T_K=peak_T,
        )
        return FieldData(
            r_mm, z_mm, T_grid,
            source_label=f"data: SYNTHETIC fallback (V1.5 tprint missing -- {e})",
            using_real=False,
        )


def mirror(r_half_mm: np.ndarray, T_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r_full = np.concatenate([-r_half_mm[:0:-1], r_half_mm])
    T_full = np.concatenate([T_grid[:, :0:-1], T_grid], axis=1)
    return r_full, T_full


def make_norm(Tmin: float, Tmax: float, gamma: float = 0.45) -> PowerNorm:
    if Tmax <= Tmin:
        Tmax = Tmin + 1.0
    return PowerNorm(gamma=gamma, vmin=Tmin, vmax=Tmax)


def zoom_window_T_range(
    T_grid: np.ndarray, r_mm: np.ndarray, z_mm: np.ndarray,
    r_lo_mm: float, r_hi_mm: float,
    z_lo_um: float, z_hi_um: float,
    mirrored: bool = False,
) -> tuple[float, float]:
    """Temperature min/max inside a zoom rectangle (for clipped color scale)."""
    z_um = z_mm * 1000.0
    j0 = int(np.searchsorted(z_um, z_lo_um, side="left"))
    j1 = int(np.searchsorted(z_um, z_hi_um, side="right"))
    j0 = max(0, j0)
    j1 = min(T_grid.shape[0], j1)

    if mirrored:
        r_full, T_full = mirror(r_mm, T_grid)
        i0 = int(np.searchsorted(r_full, r_lo_mm, side="left"))
        i1 = int(np.searchsorted(r_full, r_hi_mm, side="right"))
        i0 = max(0, i0)
        i1 = min(T_full.shape[1], i1)
        patch = T_full[j0:j1, i0:i1]
    else:
        i0 = int(np.searchsorted(r_mm, r_lo_mm, side="left"))
        i1 = int(np.searchsorted(r_mm, r_hi_mm, side="right"))
        i0 = max(0, i0)
        i1 = min(T_grid.shape[1], i1)
        patch = T_grid[j0:j1, i0:i1]

    if patch.size == 0:
        return float(T_grid.min()), float(T_grid.max())
    return float(patch.min()), float(patch.max())


# =============================================================================
# Panel drawers
# =============================================================================
def draw_half_section(
    ax, r_mm, z_mm, T_grid, norm: PowerNorm,
    *, zoom: bool,
) -> None:
    z_um = z_mm * 1000.0
    ax.pcolormesh(r_mm, z_um, T_grid, cmap="inferno", norm=norm, shading="gouraud")

    R_mm = float(r_mm[-1])
    H_um = float(z_um[-1])

    if zoom:
        ax.set_xlim(0.0, ZOOM_R_HALF_UM / 1000.0)
        ax.set_ylim(ZOOM_Z_MIN_UM, ZOOM_Z_MAX_UM)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        ax.axvline(0.0, color="cyan", linewidth=1.4, linestyle="--", alpha=0.85)
        ax.axhline(ZOOM_Z_MAX_UM, color="white", linewidth=0.9, alpha=0.85)
        ax.annotate(
            "laser spot appears\nat upper-left",
            xy=(0.0, ZOOM_Z_MAX_UM),
            xytext=(0.55 * ZOOM_R_HALF_UM, 0.5 * (ZOOM_Z_MIN_UM + ZOOM_Z_MAX_UM)),
            ha="left", va="center", fontsize=8, color="white",
            arrowprops=dict(arrowstyle="->", color="white", lw=1.3),
            bbox=dict(facecolor="black", edgecolor="white", alpha=0.70,
                      boxstyle="round,pad=0.25"),
        )
        ax.text(
            0.02, 0.04,
            "r = 0 axis at left edge",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7.5, color="white",
            bbox=dict(facecolor="black", edgecolor="white", alpha=0.65,
                      boxstyle="round,pad=0.25"),
        )
        ax.set_title(
            "Zoom -- raw half section\n"
            f"r = 0--{ZOOM_R_HALF_UM:.0f} um, z = {ZOOM_Z_MIN_UM:.0f}--{ZOOM_Z_MAX_UM:.0f} um",
            fontsize=9,
        )
    else:
        ax.set_xlim(0.0, R_mm)
        ax.set_ylim(0.0, H_um)
        ax.set_xlabel("r (mm)")
        ax.set_ylabel("z (um)")
        ax.axvline(0.0, color="cyan", linewidth=1.2, linestyle="--", alpha=0.7)
        ax.axhline(H_um, color="white", linewidth=0.7, alpha=0.7)
        ax.annotate(
            "laser spot appears\nat upper-left",
            xy=(0.0, H_um),
            xytext=(0.45 * R_mm, 0.55 * H_um),
            ha="left", va="center", fontsize=8, color="white",
            arrowprops=dict(arrowstyle="->", color="white", lw=1.3),
            bbox=dict(facecolor="black", edgecolor="white", alpha=0.65,
                      boxstyle="round,pad=0.3"),
        )
        ax.text(
            0.02, 0.04,
            "raw LS-PrePost view\nr = 0 axis at left edge",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7.5, color="white",
            bbox=dict(facecolor="black", edgecolor="white", alpha=0.55,
                      boxstyle="round,pad=0.25"),
        )
        ax.set_title(
            "Panel 1 -- raw 2-D axisymmetric half section\n"
            "(what LS-PrePost shows out of the box)",
            fontsize=9,
        )


def draw_mirrored_section(
    ax, r_mm, z_mm, T_grid, norm: PowerNorm,
    *, zoom: bool,
) -> None:
    r_full, T_full = mirror(r_mm, T_grid)
    z_um = z_mm * 1000.0
    ax.pcolormesh(r_full, z_um, T_full, cmap="inferno", norm=norm, shading="gouraud")

    R_mm = float(r_mm[-1])
    H_um = float(z_um[-1])
    r_full_um = ZOOM_R_FULL_UM / 1000.0

    if zoom:
        ax.set_xlim(-r_full_um, r_full_um)
        ax.set_ylim(ZOOM_Z_MIN_UM, ZOOM_Z_MAX_UM)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        ax.axvline(0.0, color="cyan", linewidth=1.4, linestyle="--", alpha=0.85)
        ax.axhline(ZOOM_Z_MAX_UM, color="white", linewidth=0.9, alpha=0.85)
        ax.annotate(
            "same hot zone\nafter mirror",
            xy=(0.0, ZOOM_Z_MAX_UM),
            xytext=(-0.85 * r_full_um, 0.5 * (ZOOM_Z_MIN_UM + ZOOM_Z_MAX_UM)),
            ha="left", va="center", fontsize=8, color="white",
            arrowprops=dict(arrowstyle="->", color="white", lw=1.3),
            bbox=dict(facecolor="black", edgecolor="white", alpha=0.70,
                      boxstyle="round,pad=0.25"),
        )
        ax.text(
            0.50, 0.04,
            "appears at top centre\n(this corresponds to disk centre)",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=7.5, color="white",
            bbox=dict(facecolor="black", edgecolor="cyan", alpha=0.70,
                      boxstyle="round,pad=0.25"),
        )
        ax.set_title(
            "Zoom -- mirrored full section\n"
            f"r = -{ZOOM_R_FULL_UM:.0f}--{ZOOM_R_FULL_UM:.0f} um, "
            f"z = {ZOOM_Z_MIN_UM:.0f}--{ZOOM_Z_MAX_UM:.0f} um",
            fontsize=9,
        )
    else:
        ax.set_xlim(-R_mm, R_mm)
        ax.set_ylim(0.0, H_um)
        ax.set_xlabel("r (mm)")
        ax.set_ylabel("z (um)")
        ax.axvline(0.0, color="cyan", linewidth=1.2, linestyle="--", alpha=0.7)
        ax.axhline(H_um, color="white", linewidth=0.7, alpha=0.7)
        ax.annotate(
            "same hot zone after mirror\nappears at top centre",
            xy=(0.0, H_um),
            xytext=(-0.70 * R_mm, 0.45 * H_um),
            ha="left", va="center", fontsize=8, color="white",
            arrowprops=dict(arrowstyle="->", color="white", lw=1.3),
            bbox=dict(facecolor="black", edgecolor="white", alpha=0.65,
                      boxstyle="round,pad=0.3"),
        )
        ax.text(
            0.50, 0.04,
            "this corresponds to disk centre",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=7.5, color="white",
            bbox=dict(facecolor="black", edgecolor="cyan", alpha=0.55,
                      boxstyle="round,pad=0.25"),
        )
        ax.set_title(
            "Panel 2 -- mirrored full cross-section\n"
            "(reflect r > 0 to r < 0 to match the physical disk)",
            fontsize=9,
        )


def draw_revolved_3d(
    ax, r_mm, z_mm, T_grid, norm: PowerNorm,
    *, zoom: bool,
) -> None:
    R_mm = float(r_mm[-1])
    top_T = T_grid[-1, :]

    if zoom:
        half_mm = ZOOM_DISK_HALF_UM / 1000.0
        n = 200
        xs = np.linspace(-half_mm, half_mm, n)
        ys = np.linspace(-half_mm, half_mm, n)
    else:
        n = 300
        xs = np.linspace(-R_mm, R_mm, n)
        ys = np.linspace(-R_mm, R_mm, n)

    X, Y = np.meshgrid(xs, ys, indexing="xy")
    rho = np.sqrt(X ** 2 + Y ** 2)
    T_disk = np.interp(rho, r_mm, top_T, left=top_T[0], right=top_T[-1])

    ax.pcolormesh(xs, ys, T_disk, cmap="inferno", norm=norm, shading="gouraud")

    if zoom:
        ax.set_xlim(-half_mm, half_mm)
        ax.set_ylim(-half_mm, half_mm)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v * 1000:.0f}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v * 1000:.0f}"))
        ax.set_title(
            "Zoom -- revolved top view\n"
            f"x, y = -{ZOOM_DISK_HALF_UM:.0f}--{ZOOM_DISK_HALF_UM:.0f} um",
            fontsize=9,
        )
    else:
        theta = np.linspace(0, 2 * np.pi, 360)
        ax.plot(R_mm * np.cos(theta), R_mm * np.sin(theta),
                color="white", linewidth=1.0, alpha=0.9)
        ax.set_xlim(-R_mm * 1.05, R_mm * 1.05)
        ax.set_ylim(-R_mm * 1.05, R_mm * 1.05)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(
            "Panel 3 -- 3-D revolved schematic (top view)\n"
            "Panel 1 upper-left corner -> disk centre",
            fontsize=9,
        )

    ax.set_aspect("equal")
    ax.plot([0.0], [0.0], marker="o", color="cyan", markersize=7 if zoom else 8)
    label = "laser spot at disk centre" if zoom else "laser spot\n(disk centre)"
    offset = (0.35 if zoom else 0.55) * (xs[-1] if zoom else R_mm)
    ax.annotate(
        label,
        xy=(0.0, 0.0),
        xytext=(offset, offset),
        ha="left", va="bottom", fontsize=8, color="white",
        arrowprops=dict(arrowstyle="->", color="cyan", lw=1.3),
        bbox=dict(facecolor="black", edgecolor="cyan", alpha=0.65,
                  boxstyle="round,pad=0.25"),
    )
    if not zoom:
        ax.text(
            0.02, 0.98, "revolved top view",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, color="white",
            bbox=dict(facecolor="black", edgecolor="white", alpha=0.55,
                      boxstyle="round,pad=0.25"),
        )


# =============================================================================
# Figure builders
# =============================================================================
def build_original_figure(data: FieldData, out_path: Path) -> None:
    Tmin = float(data.T_grid.min())
    Tmax = float(data.T_grid.max())
    norm = make_norm(Tmin, Tmax)

    fig = plt.figure(figsize=(16.5, 6.0))
    gs = fig.add_gridspec(1, 3, wspace=0.32)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    draw_half_section(ax1, data.r_mm, data.z_mm, data.T_grid, norm, zoom=False)
    im = ax1.collections[0] if ax1.collections else None
    draw_mirrored_section(ax2, data.r_mm, data.z_mm, data.T_grid, norm, zoom=False)
    draw_revolved_3d(ax3, data.r_mm, data.z_mm, data.T_grid, norm, zoom=False)

    cbar_ax = fig.add_axes([0.945, 0.18, 0.012, 0.68])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Temperature (K)")

    fig.suptitle(
        "V1.6 -- coordinate convention for the 2-D axisymmetric model\n"
        f"({data.source_label})",
        y=0.995, fontsize=12,
    )
    fig.subplots_adjust(left=0.05, right=0.93, top=0.88, bottom=0.10)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[OK] wrote {out_path}")


def build_zoomed_composite_figure(data: FieldData, out_path: Path) -> None:
    Tmin_full = float(data.T_grid.min())
    Tmax_full = float(data.T_grid.max())
    norm_full = make_norm(Tmin_full, Tmax_full)

    # Clipped color range from zoom windows only (hot zone visibility)
    tz_lo, tz_hi = zoom_window_T_range(
        data.T_grid, data.r_mm, data.z_mm,
        0.0, ZOOM_R_HALF_UM / 1000.0,
        ZOOM_Z_MIN_UM, ZOOM_Z_MAX_UM, mirrored=False,
    )
    tz_lo_m, tz_hi_m = zoom_window_T_range(
        data.T_grid, data.r_mm, data.z_mm,
        -ZOOM_R_FULL_UM / 1000.0, ZOOM_R_FULL_UM / 1000.0,
        ZOOM_Z_MIN_UM, ZOOM_Z_MAX_UM, mirrored=True,
    )
    # Disk zoom: use surface row within radial extent
    top_T = data.T_grid[-1, :]
    r_um = data.r_mm * 1000.0
    mask = r_um <= ZOOM_DISK_HALF_UM
    tz_disk_lo = float(top_T[mask].min()) if mask.any() else Tmin_full
    tz_disk_hi = float(top_T[mask].max()) if mask.any() else Tmax_full
    Tmin_zoom = min(tz_lo, tz_lo_m, tz_disk_lo)
    Tmax_zoom = max(tz_hi, tz_hi_m, tz_disk_hi)
    norm_zoom = make_norm(Tmin_zoom, Tmax_zoom, gamma=0.40)

    fig = plt.figure(figsize=(16.5, 10.5))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.32)

    # Row 0 -- full views
    ax00 = fig.add_subplot(gs[0, 0])
    ax01 = fig.add_subplot(gs[0, 1])
    ax02 = fig.add_subplot(gs[0, 2])
    draw_half_section(ax00, data.r_mm, data.z_mm, data.T_grid, norm_full, zoom=False)
    im_full = ax00.collections[0]
    draw_mirrored_section(ax01, data.r_mm, data.z_mm, data.T_grid, norm_full, zoom=False)
    draw_revolved_3d(ax02, data.r_mm, data.z_mm, data.T_grid, norm_full, zoom=False)
    for ax, tag in zip((ax00, ax01, ax02), ("A", "B", "C")):
        ax.text(0.02, 0.98, f"Full view {tag}", transform=ax.transAxes,
                ha="left", va="top", fontsize=8, color="white",
                bbox=dict(facecolor="0.15", edgecolor="white", alpha=0.75,
                          boxstyle="round,pad=0.2"))

    # Row 1 -- zoomed views
    ax10 = fig.add_subplot(gs[1, 0])
    ax11 = fig.add_subplot(gs[1, 1])
    ax12 = fig.add_subplot(gs[1, 2])
    draw_half_section(ax10, data.r_mm, data.z_mm, data.T_grid, norm_zoom, zoom=True)
    im_zoom = ax10.collections[0]
    draw_mirrored_section(ax11, data.r_mm, data.z_mm, data.T_grid, norm_zoom, zoom=True)
    draw_revolved_3d(ax12, data.r_mm, data.z_mm, data.T_grid, norm_zoom, zoom=True)
    for ax, tag in zip((ax10, ax11, ax12), ("A", "B", "C")):
        ax.text(0.02, 0.98, f"Zoom {tag}", transform=ax.transAxes,
                ha="left", va="top", fontsize=8, color="white",
                bbox=dict(facecolor="0.15", edgecolor="cyan", alpha=0.75,
                          boxstyle="round,pad=0.2"))

    # Two colorbars: full (left) and zoom (right)
    cbar_full_ax = fig.add_axes([0.945, 0.56, 0.012, 0.36])
    cbar_full = fig.colorbar(im_full, cax=cbar_full_ax)
    cbar_full.set_label("Temperature (K)  -- full range")

    cbar_zoom_ax = fig.add_axes([0.945, 0.12, 0.012, 0.36])
    cbar_zoom = fig.colorbar(im_zoom, cax=cbar_zoom_ax)
    cbar_zoom.set_label("Temperature (K)  -- zoomed range")

    fig.suptitle(
        "V1.6 -- axisymmetric coordinate explanation (full + zoomed views)\n"
        f"({data.source_label})\n"
        "Bottom row uses zoomed color range for visibility "
        f"({Tmin_zoom:.0f}--{Tmax_zoom:.0f} K in hot zone vs "
        f"{Tmin_full:.0f}--{Tmax_full:.0f} K full wafer)",
        y=0.995, fontsize=11,
    )
    fig.subplots_adjust(left=0.05, right=0.93, top=0.90, bottom=0.06)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[OK] wrote {out_path}")


# =============================================================================
# Main
# =============================================================================
def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--toml", type=Path,
        default=project_root / "config" / "v16_30ps.toml",
    )
    ap.add_argument(
        "--out-root", type=Path,
        default=project_root / "results" / "v16_30ps",
    )
    args = ap.parse_args()

    if not args.toml.is_file():
        print(f"[ERROR] V1.6 config not found: {args.toml}", file=sys.stderr)
        return 2
    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)

    fig_dir = args.out_root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    data = load_field_data(project_root, cfg)
    src = "real V1.5 data" if data.using_real else "synthetic fallback"
    print(f"V1.6 coordinate explanation -- {src}")

    out_orig = fig_dir / "v16_axisymmetric_coordinate_explanation.png"
    out_zoom = fig_dir / "v16_axisymmetric_coordinate_explanation_zoomed.png"
    build_original_figure(data, out_orig)
    build_zoomed_composite_figure(data, out_zoom)

    print(f"\nZoom windows:")
    print(f"  half section:     r = 0--{ZOOM_R_HALF_UM:.0f} um, "
          f"z = {ZOOM_Z_MIN_UM:.0f}--{ZOOM_Z_MAX_UM:.0f} um")
    print(f"  mirrored section: r = -{ZOOM_R_FULL_UM:.0f}--{ZOOM_R_FULL_UM:.0f} um, "
          f"z = {ZOOM_Z_MIN_UM:.0f}--{ZOOM_Z_MAX_UM:.0f} um")
    print(f"  revolved top:     x, y = -{ZOOM_DISK_HALF_UM:.0f}--{ZOOM_DISK_HALF_UM:.0f} um")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
