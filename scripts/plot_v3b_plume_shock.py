"""
plot_v3b_plume_shock.py
=======================

V3B / V3B.1 -- plot parametric plume / shock / ejecta model figures.

Run after build_v3b_plume_model.py:
    .\\.venv\\Scripts\\python.exe scripts\\plot_v3b_plume_shock.py
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import tomllib
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.patches import Rectangle

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

V3B_BANNER = "V3B.1 -- 30 ps parametric plume/shock/ejecta (driver: ep_5p0uj / V2.6)"


def read_metrics(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    bool_cols = {"melt_exists", "vapor_exists", "plume_active", "ejecta_active"}
    out: dict[str, list] = {k: [] for k in rows[0]}
    for row in rows:
        for k, v in row.items():
            if k in bool_cols:
                out[k].append(v.strip().lower() in ("yes", "true", "1"))
            else:
                out[k].append(float(v))
    return {k: np.asarray(v) for k, v in out.items()}


def metric_R0(M: dict[str, np.ndarray]) -> float:
    if "R0_um" in M:
        return float(M["R0_um"][0])
    return float(M["R_crater_um"][0])


def read_t0_ps(path: Path) -> float:
    text = path.read_text(encoding="utf-8")
    m = re.search(r"t0_ps=([0-9.eE+-]+)", text)
    if m:
        return float(m.group(1))
    return float(text.strip().splitlines()[0])


def read_v26_ep5_metrics(project_root: Path) -> dict[str, np.ndarray] | None:
    p = project_root / "results/v26_30ps_threshold_ablation/ep_5p0uj/v26_threshold_metrics.csv"
    if not p.is_file():
        return None
    return read_metrics(p)


def density_field(
    R: np.ndarray,
    Z: np.ndarray,
    sigma_r: float,
    sigma_z: float,
    z_center: float,
    n0: float,
    R0: float,
    z_surface: float,
) -> np.ndarray:
    """Semi-ellipsoid Gaussian + near-surface root term (V3B.1)."""
    if n0 <= 0:
        return np.zeros_like(R)
    rr = np.abs(R)
    zz = Z
    density = n0 * np.exp(-(rr / sigma_r) ** 2 - ((zz - z_center) / sigma_z) ** 2)
    root = np.exp(-(rr / (0.8 * R0)) ** 2) * np.exp(-(zz / 10.0) ** 2)
    density = density + 0.25 * n0 * root
    return np.where(zz >= z_surface, density, 0.0)


def semi_circle(R_val: float, z_surface: float, n: int = 120) -> tuple[np.ndarray, np.ndarray]:
    th = np.linspace(0, np.pi, n)
    return R_val * np.cos(th), z_surface + R_val * np.sin(th)


def nearest_index(M: dict[str, np.ndarray], time_ps: float) -> int:
    return int(np.argmin(np.abs(M["time_ps"] - time_ps)))


def _grid_window(cfg: dict, section: str = "domain") -> tuple[tuple[float, float], tuple[float, float], float]:
    dom = cfg[section] if section in cfg else cfg["domain"]
    r = (float(dom["r_min_um"]), float(dom["r_max_um"]))
    z = (float(dom["z_min_um"]), float(dom["z_max_um"]))
    z_surf = float(cfg["domain"]["surface_z_um"])
    return r, z, z_surf


def add_surface_bar(ax, r_win, z_win, z_surf, height: float = 8.0):
    ax.add_patch(Rectangle(
        (r_win[0], z_win[0]), r_win[1] - r_win[0], height,
        facecolor="0.25", edgecolor="none", alpha=0.9, zorder=2,
    ))
    ax.axhline(z_surf, color="white", lw=0.8, alpha=0.8)


def overlay_fronts(ax, R0, Rp, Rs, z_surf):
    if R0 > 0:
        ax.plot([-R0, R0], [z_surf, z_surf], color="#ffeb3b", lw=2.0, zorder=4)
    if Rp > 0:
        xp, yp = semi_circle(Rp, z_surf)
        ax.plot(xp, yp, color="#00bcd4", ls="--", lw=1.4, zorder=5)
    if Rs > 0:
        xs, ys = semi_circle(Rs, z_surf)
        ax.plot(xs, ys, color="#ff5252", lw=1.4, zorder=5)


def plot_plume_shock_vs_time(M, t0_ps, R0, out_path):
    t_ns = M["time_ns"]
    t0_ns = t0_ps * 1e-3
    Rp = M["R_plume_um"]
    Rs = M["R_shock_um"]

    fig, ax = plt.subplots(figsize=(11.0, 6.0))
    ax.fill_between(t_ns, Rp, Rs, where=Rs > Rp, color="#9467bd", alpha=0.25, label="shock ahead")
    ax.plot(t_ns, Rp, color="#d62728", lw=2.2, label="R_plume")
    ax.plot(t_ns, Rs, color="#9467bd", lw=2.2, label="R_shock")
    ax.axhline(R0, color="#888888", ls=":", lw=1.2, label=f"R0={R0:.2f} um")

    markers = [(t0_ps * 1e-3, "t0"), (0.1, "100 ps"), (0.5, "500 ps"), (1.0, "1 ns"), (5.0, "5 ns")]
    for t_m, lbl in markers:
        ax.axvline(t_m, color="0.55", ls=":", lw=0.8, alpha=0.7)
        ax.text(t_m, ax.get_ylim()[1] * 0.02 if ax.get_ylim()[1] > 0 else 1,
                f" {lbl}", fontsize=7, color="0.45", va="bottom")

    ax.set_xlabel("time (ns)  [1000 ps = 1 ns]")
    ax.set_ylabel("radius (um)")
    ax.set_xlim(0, float(t_ns.max()))
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.set_title(f"{V3B_BANNER}\nplume & shock front vs time  (t0 = {t0_ps:.3f} ps)")
    ax.legend(loc="upper left", fontsize=8)

    k5 = nearest_index(M, 5000)
    txt = (
        f"R0 = {R0:.2f} um\n"
        f"t0 = {t0_ps:.3f} ps\n"
        f"R_plume(5 ns) = {Rp[k5]:.1f} um\n"
        f"R_shock(5 ns) = {Rs[k5]:.1f} um\n"
        f"model = parametric proxy,\nnot CFD/SPH"
    )
    ax.text(
        0.98, 0.02, txt, transform=ax.transAxes, ha="right", va="bottom",
        fontsize=8, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_density_sequence(M, t0_ps, snapshot_times_ps, cfg, out_path):
    r_win, z_win, z_surf = _grid_window(cfg)
    R0 = metric_R0(M)
    r = np.linspace(r_win[0], r_win[1], 401)
    z = np.linspace(z_win[0], z_win[1], 301)
    R, Z = np.meshgrid(r, z, indexing="xy")

    snap_idx = [nearest_index(M, tt) for tt in snapshot_times_ps]
    fields = []
    for k in snap_idx:
        fields.append(density_field(
            R, Z,
            float(M["sigma_r_um"][k]), float(M["sigma_z_um"][k]),
            float(M["z_center_um"][k]), float(M["n0_rel"][k]), R0, z_surf,
        ))
    nmax = max(float(f.max()) for f in fields) or 1.0
    cmap = plt.get_cmap(str(cfg["render"]["colormap"]))
    norm = PowerNorm(gamma=0.45, vmin=0, vmax=nmax)

    fig, axes = plt.subplots(3, 3, figsize=(14.0, 12.5), squeeze=False)
    im = None
    for ax, k, t_tgt, field in zip(axes.ravel(), snap_idx, snapshot_times_ps, fields):
        im = ax.pcolormesh(r, z, field, cmap=cmap, norm=norm, shading="gouraud")
        add_surface_bar(ax, r_win, z_win, z_surf)
        Rp = float(M["R_plume_um"][k])
        Rs = float(M["R_shock_um"][k])
        overlay_fronts(ax, R0, Rp, Rs, z_surf)
        ax.set_xlim(r_win)
        ax.set_ylim(z_win)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        t_act = float(M["time_ps"][k])
        if t_act < t0_ps - 1e-6:
            ax.set_title(f"t={t_act:.0f} ps\nbefore plume launch", fontsize=8)
        else:
            ax.set_title(
                f"t={t_act:.0f} ps\nRplume={Rp:.1f} um, Rshock={Rs:.1f} um",
                fontsize=8,
            )

    fig.subplots_adjust(left=0.06, right=0.88, top=0.93, bottom=0.06, hspace=0.42, wspace=0.28)
    cax = fig.add_axes([0.905, 0.15, 0.018, 0.72])
    fig.colorbar(im, cax=cax, label="density proxy n(r,z)")
    fig.suptitle(
        f"{V3B_BANNER}\ndensity proxy sequence (t0={t0_ps:.2f} ps)",
        y=0.98, fontsize=11,
    )
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def plain_image_field(R, Z, field, R_plume, R_shock, z_surf, cfg):
    img = cfg["plain_imaging"]
    nmax = float(field.max()) if field.max() > 0 else 1.0
    density_norm = field / nmax
    I = img["I_bg"] * np.exp(-img["tau_eff"] * density_norm)
    I = np.where(Z >= z_surf, I, img["silicon_gray"])
    return I, density_norm


def plot_plain_imaging_sequence(M, t0_ps, snapshot_times_ps, cfg, out_path):
    r_win, z_win, z_surf = _grid_window(cfg)
    R0 = metric_R0(M)
    r = np.linspace(r_win[0], r_win[1], 401)
    z = np.linspace(z_win[0], z_win[1], 301)
    R, Z = np.meshgrid(r, z, indexing="xy")
    snap_idx = [nearest_index(M, tt) for tt in snapshot_times_ps]
    bg = float(cfg["plain_imaging"]["I_bg"])

    fig, axes = plt.subplots(3, 3, figsize=(14.0, 12.5), squeeze=False)
    for ax, k, t_tgt in zip(axes.ravel(), snap_idx, snapshot_times_ps):
        dens = density_field(
            R, Z,
            float(M["sigma_r_um"][k]), float(M["sigma_z_um"][k]),
            float(M["z_center_um"][k]), float(M["n0_rel"][k]), R0, z_surf,
        )
        I, _ = plain_image_field(
            R, Z, dens,
            float(M["R_plume_um"][k]), float(M["R_shock_um"][k]),
            z_surf, cfg,
        )
        ax.imshow(
            I, origin="lower", aspect="auto",
            extent=[r_win[0], r_win[1], z_win[0], z_win[1]],
            cmap="gray", vmin=0, vmax=bg,
        )
        ax.add_patch(Rectangle(
            (r_win[0], z_win[0]), r_win[1] - r_win[0], 8.0,
            facecolor="0.30", edgecolor="none", alpha=0.85,
        ))
        Rp, Rs = float(M["R_plume_um"][k]), float(M["R_shock_um"][k])
        if Rp > 0:
            xp, yp = semi_circle(Rp, z_surf)
            ax.plot(xp, yp, color="#00bcd4", ls="--", lw=1.2, alpha=0.9)
        if Rs > 0:
            xs, ys = semi_circle(Rs, z_surf)
            ax.plot(xs, ys, color="#ff8a80", lw=1.2, alpha=0.9)
        t_act = float(M["time_ps"][k])
        status = "" if t_act >= t0_ps else " (pre-launch)"
        ax.set_title(f"t={t_act:.0f} ps  R_plume={Rp:.0f}  R_shock={Rs:.0f}{status}", fontsize=8)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")

    tau = cfg["plain_imaging"]["tau_eff"]
    fig.suptitle(
        f"{V3B_BANNER}\nplain-imaging side view (dark = plume extinction, tau_eff={tau})",
        y=0.98, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def generate_ejecta_particles(
    n_particles: int,
    R0: float,
    ejecta_h: float,
    ejecta_r: float,
    z_surf: float,
    rng: np.random.Generator,
    brightness_scale: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if ejecta_h <= 0 or ejecta_r <= 0:
        return np.array([]), np.array([]), np.array([])
    n_eff = max(int(n_particles * brightness_scale), 0)
    if n_eff == 0:
        return np.array([]), np.array([]), np.array([])
    rr = rng.uniform(0, ejecta_r, n_eff)
    th = rng.uniform(0, 2 * np.pi, n_eff)
    x = rr * np.cos(th)
    z = z_surf + rng.uniform(0, max(ejecta_h, 1e-6), n_eff)
    return x, z, z


def ejecta_brightness_scale(t_act: float, t0_ps: float) -> float:
    if t_act < t0_ps:
        return 0.0
    t_ns = (t_act - t0_ps) * 1e-3
    rise = 1.0 - np.exp(-t_ns / 0.15)
    if t_ns <= 0.5:
        return float(rise)
    decay = np.exp(-(t_ns - 0.5) / 2.0)
    return float(max(rise * (0.35 + 0.65 * decay), 0.05))


def plot_ejecta_sequence(M, t0_ps, snapshot_times_ps, cfg, out_path):
    r_win, z_win, z_surf = _grid_window(cfg)
    ejecta_cfg = cfg["ejecta"]
    seed = int(ejecta_cfg.get("random_seed", 42))
    n_particles = int(ejecta_cfg["particle_count"])
    R0 = metric_R0(M)
    snap_idx = [nearest_index(M, tt) for tt in snapshot_times_ps]

    fig, axes = plt.subplots(3, 3, figsize=(14.0, 12.5), squeeze=False)
    for ax, k, t_tgt in zip(axes.ravel(), snap_idx, snapshot_times_ps):
        add_surface_bar(ax, r_win, z_win, z_surf, height=10.0)
        t_act = float(M["time_ps"][k])
        eh = float(M["ejecta_height_um"][k])
        er = float(M["ejecta_radius_um"][k])
        scale = ejecta_brightness_scale(t_act, t0_ps)
        if scale > 0 and eh > 0:
            sub_rng = np.random.default_rng(seed + k)
            x, z, zh = generate_ejecta_particles(
                n_particles, R0, eh, er, z_surf, sub_rng, scale,
            )
            if x.size:
                sc = ax.scatter(
                    x, z, s=6, c=zh, cmap="plasma", alpha=0.5,
                    vmin=z_surf, vmax=max(eh + z_surf, z_surf + 1),
                    edgecolors="none",
                )
        else:
            ax.text(0.5, 0.5, "no ejecta before t0", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color="0.5")
        ax.set_xlim(r_win)
        ax.set_ylim(z_win)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        ax.set_title(
            f"t={t_act:.0f} ps  h={eh:.1f} um\n"
            f"display-only ejecta proxy, not SPH particles",
            fontsize=7,
        )
    fig.suptitle(
        f"{V3B_BANNER}\nejecta particle proxy (seed={seed}, n={n_particles})",
        y=0.98, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_early_time_zoom(M, t0_ps, cfg, out_path):
    times = list(cfg["time"].get("early_zoom_times_ps", [0, 30, 50, 60, 80, 100, 150, 200, 500]))
    r_win, z_win, z_surf = _grid_window(cfg, "early_zoom")
    R0 = metric_R0(M)
    r = np.linspace(r_win[0], r_win[1], 301)
    z = np.linspace(z_win[0], z_win[1], 201)
    R, Z = np.meshgrid(r, z, indexing="xy")
    ejecta_cfg = cfg["ejecta"]
    seed = int(ejecta_cfg.get("random_seed", 42))
    n_particles = int(ejecta_cfg["particle_count"])

    fig, axes = plt.subplots(3, 3, figsize=(14.0, 12.5), squeeze=False)
    for ax, t_tgt in zip(axes.ravel(), times):
        k = nearest_index(M, t_tgt)
        field = density_field(
            R, Z,
            float(M["sigma_r_um"][k]), float(M["sigma_z_um"][k]),
            float(M["z_center_um"][k]), float(M["n0_rel"][k]), R0, z_surf,
        )
        nmax = float(field.max()) or 1.0
        ax.pcolormesh(r, z, field, cmap="inferno", shading="gouraud",
                      norm=PowerNorm(gamma=0.45, vmin=0, vmax=nmax))
        add_surface_bar(ax, r_win, z_win, z_surf, height=5.0)
        Rp = float(M["R_plume_um"][k])
        Rs = float(M["R_shock_um"][k])
        eh = float(M["ejecta_height_um"][k])
        overlay_fronts(ax, R0, Rp, Rs, z_surf)
        t_act = float(M["time_ps"][k])
        scale = ejecta_brightness_scale(t_act, t0_ps)
        if scale > 0 and eh > 0:
            sub_rng = np.random.default_rng(seed + k + 1000)
            x, pz, _ = generate_ejecta_particles(
                n_particles, R0, eh, float(M["ejecta_radius_um"][k]),
                z_surf, sub_rng, scale,
            )
            if x.size:
                ax.scatter(x, pz, s=8, c="#ff9800", alpha=0.45, edgecolors="none", zorder=6)
        ax.set_xlim(r_win)
        ax.set_ylim(z_win)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        ax.set_title(
            f"t={t_act:.1f} ps  Rplume={Rp:.1f}  Rshock={Rs:.1f}\n"
            f"ejecta h={eh:.1f} um",
            fontsize=7,
        )
    fig.suptitle(
        f"{V3B_BANNER}\nearly-time zoom 0--500 ps (r +/-80 um, z 0--120 um)",
        y=0.98, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_v4_input_fields(M, t0_ps, cfg, out_path):
    r_win, z_win, z_surf = _grid_window(cfg)
    R0 = metric_R0(M)
    r = np.linspace(r_win[0], r_win[1], 401)
    z = np.linspace(z_win[0], z_win[1], 301)
    R, Z = np.meshgrid(r, z, indexing="xy")
    tau = float(cfg["plain_imaging"]["tau_eff"])

    def frame_at(t_ps: float):
        k = nearest_index(M, t_ps)
        dens = density_field(
            R, Z,
            float(M["sigma_r_um"][k]), float(M["sigma_z_um"][k]),
            float(M["z_center_um"][k]), float(M["n0_rel"][k]), R0, z_surf,
        )
        nmax = float(dens.max()) or 1.0
        dn = dens / nmax
        return dens, dn, float(M["time_ps"][k])

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 10.0))

    dens100, dn100, t100 = frame_at(100)
    im0 = axes[0, 0].pcolormesh(r, z, dens100, cmap="inferno", shading="gouraud")
    axes[0, 0].set_title(f"density proxy @ {t100:.0f} ps")
    axes[0, 0].set_xlabel("r (um)")
    axes[0, 0].set_ylabel("z (um)")
    fig.colorbar(im0, ax=axes[0, 0], shrink=0.85, label="n proxy")

    dens500, dn500, t500 = frame_at(500)
    im1 = axes[0, 1].pcolormesh(r, z, dens500, cmap="inferno", shading="gouraud")
    axes[0, 1].set_title(f"density proxy @ {t500:.0f} ps")
    axes[0, 1].set_xlabel("r (um)")
    axes[0, 1].set_ylabel("z (um)")
    fig.colorbar(im1, ax=axes[0, 1], shrink=0.85, label="n proxy")

    B = np.exp(-tau * dn500)
    im2 = axes[1, 0].pcolormesh(r, z, B, cmap="gray_r", shading="gouraud", vmin=0, vmax=1)
    axes[1, 0].set_title(f"modulation loss B=exp(-tau*n) @ {t500:.0f} ps  (tau={tau})")
    axes[1, 0].set_xlabel("r (um)")
    axes[1, 0].set_ylabel("z (um)")
    fig.colorbar(im2, ax=axes[1, 0], shrink=0.85, label="B proxy")

    dr, dz = np.gradient(dn500, r[1] - r[0], z[1] - z[0])
    phi_proxy = np.cumsum(dz, axis=0)
    phi_proxy = phi_proxy / (float(phi_proxy.max()) or 1.0)
    im3 = axes[1, 1].pcolormesh(r, z, phi_proxy, cmap="RdBu_r", shading="gouraud")
    axes[1, 1].set_title(f"phase distortion proxy (grad integral) @ {t500:.0f} ps")
    axes[1, 1].set_xlabel("r (um)")
    axes[1, 1].set_ylabel("z (um)")
    fig.colorbar(im3, ax=axes[1, 1], shrink=0.85, label="phi proxy")

    for ax in axes.ravel():
        ax.set_xlim(r_win)
        ax.set_ylim(z_win)

    fig.suptitle(
        f"{V3B_BANNER}\nV4 input field check (proxy only, not physical phase)",
        y=0.98, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def plot_driver_summary(M, t0_ps, cfg, v26_m, out_path):
    R0 = metric_R0(M)
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 10.0))

    ax = axes[0, 0]
    if v26_m is not None:
        t = v26_m["time_ps"]
        mask = t <= 120
        ax.plot(t[mask], v26_m["crater_depth_um"][mask], "r-", lw=1.8, label="crater depth")
        ax.plot(t[mask], v26_m["crater_radius_um"][mask], "b--", lw=1.2, label="crater radius")
        ax.axvline(30.029, color="#2ca02c", ls=":", lw=1.0)
        ax.axvline(50.929, color="#9467bd", ls=":", lw=1.0)
    ax.set_xlim(0, 120)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("um")
    ax.set_title("V2.6 driver (ep_5p0uj, 0--120 ps)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    mask = M["time_ns"] <= 5.0
    ax.fill_between(
        M["time_ns"][mask], M["R_plume_um"][mask], M["R_shock_um"][mask],
        where=M["R_shock_um"][mask] > M["R_plume_um"][mask],
        color="#9467bd", alpha=0.2,
    )
    ax.plot(M["time_ns"][mask], M["R_plume_um"][mask], "r-", lw=1.8, label="R_plume")
    ax.plot(M["time_ns"][mask], M["R_shock_um"][mask], color="#9467bd", lw=1.8, label="R_shock")
    ax.axvline(t0_ps * 1e-3, color="0.35", ls="--", lw=1.0)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("radius (um)")
    ax.set_title("V3B.1 fronts (0--5 ns)")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    r_win, z_win, z_surf = _grid_window(cfg)
    r = np.linspace(r_win[0], r_win[1], 301)
    z = np.linspace(z_win[0], z_win[1], 201)
    R, Z = np.meshgrid(r, z, indexing="xy")
    k0 = nearest_index(M, 500)
    field = density_field(
        R, Z,
        float(M["sigma_r_um"][k0]), float(M["sigma_z_um"][k0]),
        float(M["z_center_um"][k0]), float(M["n0_rel"][k0]), R0, z_surf,
    )
    ax = axes[1, 0]
    im = ax.pcolormesh(r, z, field, cmap="inferno", shading="gouraud")
    overlay_fronts(ax, R0, float(M["R_plume_um"][k0]), float(M["R_shock_um"][k0]), z_surf)
    ax.set_xlim(r_win)
    ax.set_ylim(0, 250)
    ax.set_title(f"density @ {float(M['time_ps'][k0]):.0f} ps")
    ax.set_xlabel("r (um)")
    ax.set_ylabel("z (um)")
    fig.colorbar(im, ax=ax, shrink=0.85, label="n proxy")

    ax = axes[1, 1]
    k_rep = nearest_index(M, 200)
    eh = float(M["ejecta_height_um"][k_rep])
    er = float(M["ejecta_radius_um"][k_rep])
    rng = np.random.default_rng(int(cfg["ejecta"]["random_seed"]) + k_rep)
    x, pz, _ = generate_ejecta_particles(
        int(cfg["ejecta"]["particle_count"]), R0, eh, er, z_surf, rng, 1.0,
    )
    ax.add_patch(Rectangle((r_win[0], 0), r_win[1] - r_win[0], 10, facecolor="0.3"))
    if x.size:
        ax.scatter(x, pz, s=5, c=pz, cmap="plasma", alpha=0.55, edgecolors="none")
    ax.set_xlim(-100, 100)
    ax.set_ylim(0, 150)
    ax.set_title(f"ejecta proxy @ {float(M['time_ps'][k_rep]):.0f} ps")
    ax.set_xlabel("r (um)")
    ax.set_ylabel("z (um)")

    fig.suptitle(f"{V3B_BANNER}\nV3B.1 driver summary (V2.6 -> V3B handoff)", y=0.98, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--toml", type=Path, default=project_root / "config/plume_model_v3b_30ps.toml")
    ap.add_argument("--v3b-root", type=Path, default=project_root / "results/v3b_30ps_plume")
    args = ap.parse_args()

    csv_path = args.v3b_root / "v3b_plume_shock_metrics.csv"
    t0_path = args.v3b_root / "v3b_t0.txt"
    if not csv_path.is_file():
        print(f"[ERROR] run build first: {csv_path}", file=sys.stderr)
        return 2

    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)

    M = read_metrics(csv_path)
    t0_ps = read_t0_ps(t0_path)
    R0 = metric_R0(M)
    snapshots = list(cfg["time"]["snapshot_times_ps"])
    fig_dir = args.v3b_root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_plume_shock_vs_time(M, t0_ps, R0, fig_dir / "v3b_plume_shock_front_vs_time.png")
    plot_density_sequence(M, t0_ps, snapshots, cfg, fig_dir / "v3b_density_proxy_sequence.png")
    plot_plain_imaging_sequence(M, t0_ps, snapshots, cfg, fig_dir / "v3b_plain_imaging_sequence.png")
    plot_ejecta_sequence(M, t0_ps, snapshots, cfg, fig_dir / "v3b_ejecta_particle_proxy_sequence.png")
    plot_early_time_zoom(M, t0_ps, cfg, fig_dir / "v3b_early_time_zoom_sequence.png")
    plot_v4_input_fields(M, t0_ps, cfg, fig_dir / "v3b_v4_input_fields.png")
    plot_driver_summary(M, t0_ps, cfg, read_v26_ep5_metrics(project_root),
                        fig_dir / "v3b_driver_summary.png")

    print(f"\n[DONE] V3B.1 figures in {fig_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
