"""
plot_v43_paper_style_evolution.py
=================================

Create paper-style proxy evolution figures for V4B/V4C ablation results.

This script does not parse d3plot binary fields.  It builds a clear proxy
visualization from V4C metadata and LS-PrePost peak readings recorded in
config/v43_paper_style_evolution.toml.

Run from project root:
    python scripts\\plot_v43_paper_style_evolution.py
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import tomllib
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm, Normalize  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


def load_toml(path: Path) -> dict:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def observed_peaks(cfg: dict) -> dict[str, dict]:
    out = {}
    for row in cfg.get("observed_peaks", []):
        out[str(row["case"])] = {
            "recoil_velocity_m_s": float(row["recoil_velocity_m_s"]),
            "max_displacement_mm": float(row["max_displacement_mm"]),
            "max_pressure_GPa": float(row["max_pressure_GPa"]),
            "min_pressure_GPa": float(row["min_pressure_GPa"]),
        }
    return out


def crater_profile(x_um: np.ndarray, radius_um: float, depth_um: float) -> np.ndarray:
    z = np.zeros_like(x_um)
    inside = np.abs(x_um) <= radius_um
    z[inside] = -depth_um * (1.0 - (x_um[inside] / max(radius_um, 1e-9)) ** 2)
    return z


def substrate_field(x_um: np.ndarray, z_um: np.ndarray, t_ps: float,
                    radius_um: float, pressure_peak: float) -> np.ndarray:
    X, Z = np.meshgrid(x_um, z_um)
    # A simple wave-like proxy launched from the crater rim into the solid.
    r = np.sqrt((np.abs(X) - 0.35 * radius_um) ** 2 + (Z + 0.05) ** 2)
    front = 0.18 * t_ps
    compress = np.exp(-((r - front) / 7.5) ** 2)
    tensile = -0.55 * np.exp(-((r - front * 1.35) / 10.0) ** 2)
    near = 0.8 * np.exp(-(X / max(radius_um * 0.9, 1e-9)) ** 2) * np.exp((Z) / 1.2)
    return pressure_peak * (compress + tensile + 0.35 * near)


def make_particles(rng: np.random.Generator, n: int, radius_um: float) -> dict[str, np.ndarray]:
    u = rng.uniform(-1.0, 1.0, n)
    x0 = radius_um * np.sign(u) * np.sqrt(np.abs(u))
    angle = rng.normal(0.0, 0.33, n)
    speed = rng.lognormal(mean=0.0, sigma=0.35, size=n)
    return {"x0": x0, "angle": angle, "speed": speed}


def particle_state(base: dict[str, np.ndarray], t_ps: float, plume_velocity_m_s: float,
                   spread: float, temp_peak: float, temp_floor: float,
                   cooling_ps: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # 1 m/s = 1e-6 um/ps.
    v_um_ps = plume_velocity_m_s * 1e-6
    rise = v_um_ps * t_ps * base["speed"]
    x = base["x0"] + spread * rise * np.sin(base["angle"])
    z = rise * np.cos(base["angle"]) + 0.18 * np.sqrt(np.maximum(t_ps, 0.0))
    temp = temp_floor + (temp_peak - temp_floor) * np.exp(-t_ps / cooling_ps) * (
        0.35 + 0.65 * np.clip(base["speed"] / 2.0, 0.0, 1.0)
    )
    return x, z, temp


def plot_snapshot_sequence(cfg: dict, out_path: Path, frame_csv: Path) -> None:
    geom = cfg["geometry"]
    vis = cfg["visualization"]
    plume = cfg["plume_proxy"]
    peaks = observed_peaks(cfg)
    case = str(vis["case"])
    peak = peaks[case]
    times = [float(t) for t in vis["times_ps"]]
    radius_um = 17.388935
    depth_um = float(geom["crater_depth_um"])

    rng = np.random.default_rng(int(vis["random_seed"]))
    particles = make_particles(rng, int(plume["particle_count"]), radius_um)

    x = np.linspace(-42.0, 42.0, 360)
    z = np.linspace(-5.0, 0.05, 140)
    surface = crater_profile(x, radius_um, depth_um)

    fig, axes = plt.subplots(2, len(times), figsize=(2.35 * len(times), 5.2), sharex=True, sharey=False)
    cmap = "turbo"
    norm = LogNorm(vmin=float(plume["temperature_floor_K"]), vmax=float(plume["temperature_peak_K"]))
    frame_rows = []

    for col, t_ps in enumerate(times):
        ax = axes[0, col]
        pressure = substrate_field(x, z, t_ps, radius_um, peak["max_pressure_GPa"])
        ax.contourf(x, z, pressure, levels=18, cmap="coolwarm",
                    norm=Normalize(vmin=peak["min_pressure_GPa"], vmax=peak["max_pressure_GPa"]))
        ax.fill_between(x, -5.0, surface, color="#7367c9", alpha=0.38, lw=0)
        ax.plot(x, surface, color="white", lw=1.0)
        xp, zp, tp = particle_state(
            particles,
            t_ps,
            float(plume["plume_velocity_m_s"]),
            float(plume["radial_spread_fraction"]),
            float(plume["temperature_peak_K"]),
            float(plume["temperature_floor_K"]),
            float(plume["cooling_time_ps"]),
        )
        visible = (zp >= -0.05) & (zp <= 70.0) & (np.abs(xp) <= 45.0)
        ax.scatter(xp[visible], zp[visible], c=tp[visible], s=2.2, cmap=cmap, norm=norm, alpha=0.82, linewidths=0)
        ax.set_title(f"{t_ps:g} ps", fontsize=9)
        ax.set_xlim(-42, 42)
        ax.set_ylim(-4.0, 68.0)
        ax.set_xticks([])
        ax.set_yticks([])
        if col == 0:
            ax.set_ylabel("plume + pressure", fontsize=9)

        ax2 = axes[1, col]
        zoom_z = np.linspace(-0.65, 0.08, 120)
        pressure_zoom = substrate_field(x, zoom_z, t_ps, radius_um, peak["max_pressure_GPa"])
        ax2.contourf(x, zoom_z, pressure_zoom, levels=18, cmap="coolwarm",
                     norm=Normalize(vmin=peak["min_pressure_GPa"], vmax=peak["max_pressure_GPa"]))
        ax2.fill_between(x, -0.65, surface, color="#7967c8", alpha=0.55, lw=0)
        ax2.plot(x, surface, color="white", lw=1.1)
        ax2.set_xlim(-42, 42)
        ax2.set_ylim(-0.55, 0.12)
        ax2.set_xticks([])
        ax2.set_yticks([])
        if col == 0:
            ax2.set_ylabel("crater zoom", fontsize=9)

        frame_rows.append({
            "time_ps": f"{t_ps:g}",
            "plume_height_um_proxy": f"{np.percentile(zp[visible], 95) if np.any(visible) else 0.0:.4f}",
            "particle_temperature_max_K_proxy": f"{np.max(tp[visible]) if np.any(visible) else 0.0:.2f}",
            "pressure_max_GPa_reference": f"{peak['max_pressure_GPa']:.6g}",
            "pressure_min_GPa_reference": f"{peak['min_pressure_GPa']:.6g}",
            "max_displacement_mm_reference": f"{peak['max_displacement_mm']:.6g}",
        })

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.025, pad=0.012)
    cbar.set_label("proxy particle T (K)")
    fig.suptitle(
        "V4D proxy evolution: V4C ep_5p0uj recoil 5000 m/s\n"
        "pre-deleted crater + recoil-driven pressure/plume proxy",
        fontsize=12,
    )
    fig.subplots_adjust(left=0.035, right=0.90, top=0.84, bottom=0.06, wspace=0.03, hspace=0.08)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)

    with frame_csv.open("w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(frame_rows[0].keys()))
        wr.writeheader()
        wr.writerows(frame_rows)


def plot_scaling_summary(cfg: dict, summary_rows: list[dict[str, str]], out_path: Path) -> None:
    peaks = observed_peaks(cfg)
    data = [peaks[k] for k in sorted(peaks, key=lambda s: abs(peaks[s]["recoil_velocity_m_s"]))]
    recoil = np.array([abs(d["recoil_velocity_m_s"]) for d in data], dtype=float)
    disp_nm = np.array([d["max_displacement_mm"] * 1e6 for d in data], dtype=float)
    pmax = np.array([d["max_pressure_GPa"] for d in data], dtype=float)
    pmin = np.array([abs(d["min_pressure_GPa"]) for d in data], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    axes[0].plot(recoil, disp_nm, marker="o", color="#d62728", lw=2)
    axes[0].set_xlabel("|recoil velocity| (m/s)")
    axes[0].set_ylabel("max resultant displacement (nm)")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(recoil, pmax, marker="o", color="#ff7f0e", lw=2, label="compression")
    axes[1].plot(recoil, pmin, marker="s", color="#1f77b4", lw=2, label="tension magnitude")
    axes[1].set_xlabel("|recoil velocity| (m/s)")
    axes[1].set_ylabel("pressure magnitude (GPa)")
    axes[1].legend(frameon=False)
    axes[1].grid(True, alpha=0.3)

    labels = []
    status = []
    for row in summary_rows:
        if "recoil_" not in row.get("case", ""):
            continue
        labels.append(row["case"].replace("ep_5p0uj_", "").replace("_", "\n"))
        status.append(1 if row.get("normal_termination") == "yes" else 0)
    axes[2].bar(np.arange(len(labels)), status, color=["#2ca02c" if s else "#d62728" for s in status])
    axes[2].set_xticks(np.arange(len(labels)), labels, fontsize=8)
    axes[2].set_ylim(0, 1.2)
    axes[2].set_yticks([0, 1], ["not run/fail", "normal"])
    axes[2].set_title("run status")

    fig.suptitle("V4C recoil sweep summary from LS-PrePost peak readings", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--toml", type=Path, default=project_root / "config" / "v43_paper_style_evolution.toml")
    args = ap.parse_args()
    cfg = load_toml(args.toml)
    out_dir = project_root / cfg["visualization"]["out_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = read_csv_rows(project_root / cfg["source"]["v42_summary_csv"])
    snapshot = out_dir / "v43_snapshot_sequence.png"
    scaling = out_dir / "v43_recoil_scaling_summary.png"
    frame_csv = out_dir / "v43_proxy_frame_data.csv"
    plot_snapshot_sequence(cfg, snapshot, frame_csv)
    plot_scaling_summary(cfg, summary, scaling)
    print(f"[OK] {snapshot}")
    print(f"[OK] {scaling}")
    print(f"[OK] {frame_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
