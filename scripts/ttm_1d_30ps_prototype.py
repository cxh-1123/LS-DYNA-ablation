"""
ttm_1d_30ps_prototype.py
========================

V1.6 step 2 -- 1D two-temperature model (TTM) prototype for a 30 ps
silicon ablation pulse.

This is a *qualitative scale-finder*: the electron heat capacity
constant, the electron thermal conductivity, the electron-phonon
coupling factor, and the absorption depths are placeholder values for
silicon in the UV/visible.  See `docs/README_v16_30ps.md` for the
caveat list.  The script's purpose is to demonstrate that

    * the electron temperature Te(z, t) spikes during the 30 ps pulse,
    * the lattice temperature Tl(z, t) lags by an electron-phonon
      coupling timescale tau_ep ~ 10 - 100 ps,
    * Te can transiently exceed T_v while Tl stays below T_m,
    * the TTM equilibration timescale and the diffusion length are
      both *much* smaller than the V1 mesh cell.

It does NOT predict ablation thresholds quantitatively.

Numerics:
    Explicit forward-Euler finite-difference in z, with Neumann
    (insulated) boundary conditions at both ends.  The default
    dt = 2 fs satisfies the diffusion CFL for the Te equation at room
    temperature (where Ce is smallest); a soft warning is printed if
    instability is detected.

Inputs:
    - config/v16_30ps.toml

Outputs:
    - results/v16_30ps/ttm_surface_temperature.csv
    - results/v16_30ps/ttm_depth_profiles.csv
    - results/v16_30ps/figures/v16_ttm_surface_temperature.png
    - results/v16_30ps/figures/v16_ttm_depth_profiles.png
    - results/v16_30ps/figures/v16_ttm_threshold_map.png

Run from project root:
    .\\.venv\\Scripts\\python.exe scripts\\ttm_1d_30ps_prototype.py
"""

from __future__ import annotations

import argparse
import csv
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "savefig.dpi": 300,
})


# =============================================================================
# Config dataclass
# =============================================================================
@dataclass
class TTMConfig:
    # Material
    T0: float
    Tm: float
    Tv: float
    rho: float
    cp: float
    kl: float
    gamma_e: float
    ke: float
    G: float
    Te_clip: float
    Tl_clip: float
    # Laser
    tau_s: float
    w0_m: float
    A: float
    Ep_J_list: list[float]
    Ep_uJ_list: list[float]
    delta_abs_m: float
    wavelength_nm: float
    # Solver
    z_max_m: float
    dz_m: float
    t_start_s: float
    t_end_s: float
    dt_s: float
    snapshot_times_s: list[float]
    save_every_n: int


def load_config(toml_path: Path) -> TTMConfig:
    with toml_path.open("rb") as fh:
        cfg = tomllib.load(fh)
    mat = cfg["material"]
    ttm = mat["ttm_prototype"]
    tt  = cfg["ttm_solver"]
    sc  = cfg["ttm_scan"]
    abs_table = cfg["laser"]["absorption_depth_nm"]
    wl_nm = float(sc["wavelength_nm"])
    delta_abs_nm = float(abs_table[str(int(wl_nm))])
    return TTMConfig(
        T0           = float(mat["T0_K"]),
        Tm           = float(mat["T_melt_K"]),
        Tv           = float(mat["T_vap_K"]),
        rho          = float(mat["rho_kg_m3"]),
        cp           = float(mat["cp_J_kgK"]),
        kl           = float(mat["k_W_mK"]),
        gamma_e      = float(ttm["gamma_e_J_m3_K2"]),
        ke           = float(ttm["ke_W_mK"]),
        G            = float(ttm["G_W_m3_K"]),
        Te_clip      = float(ttm["Te_max_clip_K"]),
        Tl_clip      = float(ttm["Tl_max_clip_K"]),
        tau_s        = float(sc["pulse_width_ps"]) * 1e-12,
        w0_m         = float(sc["spot_radius_um"]) * 1e-6,
        A            = float(sc["absorption_A"]),
        Ep_uJ_list   = list(sc["pulse_energy_uJ"]),
        Ep_J_list    = [e * 1e-6 for e in sc["pulse_energy_uJ"]],
        delta_abs_m  = delta_abs_nm * 1e-9,
        wavelength_nm= wl_nm,
        z_max_m      = float(tt["z_max_nm"]) * 1e-9,
        dz_m         = float(tt["dz_nm"]) * 1e-9,
        t_start_s    = float(tt["t_start_ps"]) * 1e-12,
        t_end_s      = float(tt["t_end_ps"]) * 1e-12,
        dt_s         = float(tt["dt_target_fs"]) * 1e-15,
        snapshot_times_s = [float(t) * 1e-12 for t in tt["snapshot_times_ps"]],
        save_every_n = int(tt["save_every_n_steps"]),
    )


# =============================================================================
# TTM solver  (explicit FD)
# =============================================================================
def run_ttm_1d(cfg: TTMConfig, Ep_J: float) -> dict:
    nz = int(round(cfg.z_max_m / cfg.dz_m)) + 1
    z = np.linspace(0.0, cfg.z_max_m, nz)

    Te = np.full(nz, cfg.T0)
    Tl = np.full(nz, cfg.T0)
    cl = cfg.rho * cfg.cp     # J/(m^3 K)

    F_peak = 2.0 * Ep_J / (np.pi * cfg.w0_m ** 2)     # J/m^2

    # CFL check at the start (worst case at T = T0):
    Ce_min = max(cfg.gamma_e * cfg.T0, 1.0)
    alpha_e_max = cfg.ke / Ce_min
    dt_cfl = cfg.dz_m ** 2 / (2.0 * alpha_e_max)
    warn = []
    dt = cfg.dt_s
    if dt > 0.95 * dt_cfl:
        dt = 0.5 * dt_cfl
        warn.append(f"reduced dt from {cfg.dt_s*1e15:.2f} fs to "
                    f"{dt*1e15:.2f} fs to satisfy CFL")

    n_steps = int(np.ceil((cfg.t_end_s - cfg.t_start_s) / dt))
    dt_actual = (cfg.t_end_s - cfg.t_start_s) / n_steps

    # Pre-compute coefficients
    inv_dz2 = 1.0 / cfg.dz_m ** 2

    surface_t   : list[float] = []
    surface_Te  : list[float] = []
    surface_Tl  : list[float] = []
    Te_global   : list[float] = []
    Tl_global   : list[float] = []

    snap_times = sorted(cfg.snapshot_times_s)
    snap_next  = 0
    snap_profiles: dict[float, tuple[np.ndarray, np.ndarray]] = {}

    t = cfg.t_start_s
    runaway = False
    for step in range(n_steps):
        # Gaussian pulse temporal profile (centre at t=0)
        q_surface = (
            F_peak / (np.sqrt(2 * np.pi) * cfg.tau_s)
            * np.exp(-(t ** 2) / (2.0 * cfg.tau_s ** 2))
        )                                                   # W/m^2
        S = cfg.A * q_surface / cfg.delta_abs_m * np.exp(-z / cfg.delta_abs_m)   # W/m^3

        Ce = cfg.gamma_e * Te                               # J/(m^3 K)

        d2Te = np.empty(nz)
        d2Te[1:-1] = (Te[2:] - 2.0 * Te[1:-1] + Te[:-2]) * inv_dz2
        d2Te[0]    = 2.0 * (Te[1]  - Te[0])  * inv_dz2
        d2Te[-1]   = 2.0 * (Te[-2] - Te[-1]) * inv_dz2

        d2Tl = np.empty(nz)
        d2Tl[1:-1] = (Tl[2:] - 2.0 * Tl[1:-1] + Tl[:-2]) * inv_dz2
        d2Tl[0]    = 2.0 * (Tl[1]  - Tl[0])  * inv_dz2
        d2Tl[-1]   = 2.0 * (Tl[-2] - Tl[-1]) * inv_dz2

        coupling = cfg.G * (Te - Tl)

        Te = Te + dt_actual * (cfg.ke * d2Te - coupling + S) / Ce
        Tl = Tl + dt_actual * (cfg.kl * d2Tl + coupling)    / cl

        # Numerical safety
        if not np.isfinite(Te).all() or not np.isfinite(Tl).all():
            runaway = True
            warn.append(f"non-finite T detected at step {step}, t={t*1e12:.2f} ps. "
                        "Clipping and continuing.")
            Te = np.nan_to_num(Te, nan=cfg.Te_clip, posinf=cfg.Te_clip, neginf=cfg.T0)
            Tl = np.nan_to_num(Tl, nan=cfg.Tl_clip, posinf=cfg.Tl_clip, neginf=cfg.T0)
        if Te.max() > cfg.Te_clip:
            if not runaway:
                warn.append(f"Te exceeded {cfg.Te_clip:.0e} K at step {step}; "
                            "clipping (prototype only).")
                runaway = True
            Te = np.clip(Te, cfg.T0, cfg.Te_clip)
        if Tl.max() > cfg.Tl_clip:
            warn.append(f"Tl exceeded {cfg.Tl_clip:.0e} K at step {step}; clipping.")
            Tl = np.clip(Tl, cfg.T0, cfg.Tl_clip)
        Te = np.maximum(Te, cfg.T0)
        Tl = np.maximum(Tl, cfg.T0)

        # Sample history
        if step % cfg.save_every_n == 0:
            surface_t.append(t)
            surface_Te.append(float(Te[0]))
            surface_Tl.append(float(Tl[0]))
            Te_global.append(float(Te.max()))
            Tl_global.append(float(Tl.max()))

        while snap_next < len(snap_times) and (t + dt_actual) >= snap_times[snap_next]:
            snap_profiles[snap_times[snap_next]] = (Te.copy(), Tl.copy())
            snap_next += 1

        t += dt_actual

    return {
        "z_m":            z,
        "t_s":            np.asarray(surface_t),
        "Te_surface_K":   np.asarray(surface_Te),
        "Tl_surface_K":   np.asarray(surface_Tl),
        "Te_max_K":       np.asarray(Te_global),
        "Tl_max_K":       np.asarray(Tl_global),
        "snapshot_profiles": snap_profiles,
        "warnings":       warn,
        "dt_actual_s":    dt_actual,
        "n_steps":        n_steps,
    }


# =============================================================================
# CSV writers
# =============================================================================
def write_surface_csv(per_case: dict[float, dict], path: Path) -> None:
    rows: list[dict] = []
    for Ep_uJ, res in per_case.items():
        for k in range(res["t_s"].shape[0]):
            rows.append({
                "Ep_uJ":        f"{Ep_uJ:g}",
                "t_ps":         f"{res['t_s'][k]*1e12:.4f}",
                "Te_surface_K": f"{res['Te_surface_K'][k]:.3f}",
                "Tl_surface_K": f"{res['Tl_surface_K'][k]:.3f}",
                "Te_max_K":     f"{res['Te_max_K'][k]:.3f}",
                "Tl_max_K":     f"{res['Tl_max_K'][k]:.3f}",
            })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


def write_depth_csv(per_case: dict[float, dict], z_m: np.ndarray, path: Path) -> None:
    rows: list[dict] = []
    for Ep_uJ, res in per_case.items():
        for t_s, (Te, Tl) in sorted(res["snapshot_profiles"].items()):
            for j in range(z_m.shape[0]):
                rows.append({
                    "Ep_uJ":  f"{Ep_uJ:g}",
                    "t_ps":   f"{t_s*1e12:.3f}",
                    "z_nm":   f"{z_m[j]*1e9:.2f}",
                    "Te_K":   f"{Te[j]:.3f}",
                    "Tl_K":   f"{Tl[j]:.3f}",
                })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


# =============================================================================
# Plotting
# =============================================================================
def plot_surface_history(per_case: dict[float, dict], cfg: TTMConfig,
                         out_path: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.5))

    cmap = plt.cm.viridis(np.linspace(0.15, 0.95, len(per_case)))
    for color, (Ep_uJ, res) in zip(cmap, per_case.items()):
        t_ps = res["t_s"] * 1e12
        ax1.plot(t_ps, res["Te_surface_K"], color=color, linewidth=1.5,
                 label=f"{Ep_uJ:g} uJ")
        ax2.plot(t_ps, res["Tl_surface_K"], color=color, linewidth=1.5,
                 label=f"{Ep_uJ:g} uJ")

    for ax, name in ((ax1, "Te"), (ax2, "Tl")):
        ax.axhline(cfg.Tm, color="#00aaff", linestyle="--", linewidth=1.0,
                   label=f"T_m = {cfg.Tm:.0f} K")
        ax.axhline(cfg.Tv, color="#ff5252", linestyle="--", linewidth=1.0,
                   label=f"T_v = {cfg.Tv:.0f} K")
        ax.axvspan(-cfg.tau_s * 1e12, cfg.tau_s * 1e12,
                   color="0.85", alpha=0.4, label="pulse FWHM ~ tau")
        ax.set_xlabel("time (ps)")
        ax.set_ylabel(f"{name}(z=0, t)   (K)")
        ax.grid(True, alpha=0.3)
        ax.set_yscale("log")
        ax.set_ylim(bottom=cfg.T0 * 0.9)
        ax.legend(fontsize=7, ncol=2, loc="upper right")
        ax.set_title(f"{name} at the surface  z = 0")

    fig.suptitle(
        f"V1.6 -- 1D TTM prototype, tau = {cfg.tau_s*1e12:.0f} ps, "
        f"w0 = {cfg.w0_m*1e6:.0f} um, A = {cfg.A:.2f}, "
        f"lambda = {cfg.wavelength_nm:.0f} nm   (PROTOTYPE PARAMETERS)",
        y=1.0, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[OK] wrote {out_path}")


def plot_depth_profiles(per_case: dict[float, dict], z_m: np.ndarray,
                        cfg: TTMConfig, out_path: Path) -> None:
    # Pick the two most informative Ep values: largest sub-threshold and
    # smallest super-threshold (by Te_max).  Fall back to the first/last.
    sorted_Ep = sorted(per_case.keys())
    pick = [sorted_Ep[max(0, len(sorted_Ep) // 2 - 1)],
            sorted_Ep[-1]]
    pick = list(dict.fromkeys(pick))  # dedupe preserving order

    n_panels = len(pick)
    fig, axes = plt.subplots(2, n_panels, figsize=(7.0 * n_panels, 8.0),
                             squeeze=False)
    for col, Ep_uJ in enumerate(pick):
        res = per_case[Ep_uJ]
        snap_items = sorted(res["snapshot_profiles"].items())
        colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(snap_items)))

        ax_e = axes[0, col]
        ax_l = axes[1, col]
        for color, (t_s, (Te_z, Tl_z)) in zip(colors, snap_items):
            ax_e.plot(z_m * 1e9, Te_z, color=color, linewidth=1.5,
                      label=f"t = {t_s*1e12:.0f} ps")
            ax_l.plot(z_m * 1e9, Tl_z, color=color, linewidth=1.5,
                      label=f"t = {t_s*1e12:.0f} ps")

        for ax, name in ((ax_e, "Te"), (ax_l, "Tl")):
            ax.axhline(cfg.Tm, color="#00aaff", linestyle="--", linewidth=1.0)
            ax.axhline(cfg.Tv, color="#ff5252", linestyle="--", linewidth=1.0)
            ax.set_xlabel("z (nm)  -- 0 = surface")
            ax.set_ylabel(f"{name}(z, t)   (K)")
            ax.set_yscale("log")
            ax.set_ylim(bottom=cfg.T0 * 0.9)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=7, loc="upper right")

        axes[0, col].set_title(
            f"Ep = {Ep_uJ:g} uJ   Te depth profiles", fontsize=11,
        )
        axes[1, col].set_title(
            f"Ep = {Ep_uJ:g} uJ   Tl depth profiles", fontsize=11,
        )

    fig.suptitle(
        "V1.6 -- 1D TTM prototype:  Te(z, t) (top row) vs Tl(z, t) (bottom row)\n"
        "(PROTOTYPE PARAMETERS; cyan dashed = T_m, red dashed = T_v)",
        y=1.0, fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[OK] wrote {out_path}")


def plot_threshold_map(per_case: dict[float, dict], cfg: TTMConfig,
                       out_path: Path) -> None:
    Ep_uJ = np.array(list(per_case.keys()))
    Te_pk = np.array([per_case[e]["Te_max_K"].max() for e in Ep_uJ])
    Tl_pk = np.array([per_case[e]["Tl_max_K"].max() for e in Ep_uJ])

    fig, ax = plt.subplots(figsize=(9.0, 5.5))
    ax.plot(Ep_uJ, Te_pk, marker="o", linewidth=1.8, color="#1f77b4",
            label="Te peak (electron)")
    ax.plot(Ep_uJ, Tl_pk, marker="s", linewidth=1.8, color="#d62728",
            label="Tl peak (lattice)")

    ax.axhline(cfg.Tm, color="#00aaff", linestyle="--", linewidth=1.2,
               label=f"T_m = {cfg.Tm:.0f} K")
    ax.axhline(cfg.Tv, color="#ff5252", linestyle="--", linewidth=1.2,
               label=f"T_v = {cfg.Tv:.0f} K")

    ax.set_xlabel("pulse energy Ep (uJ)")
    ax.set_ylabel("peak temperature anywhere in domain (K)")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_title(
        f"V1.6 -- TTM peak Te / Tl vs Ep  (tau = {cfg.tau_s*1e12:.0f} ps, "
        f"w0 = {cfg.w0_m*1e6:.0f} um, A = {cfg.A:.2f}, "
        f"lambda = {cfg.wavelength_nm:.0f} nm)\n"
        "PROTOTYPE PARAMETERS -- qualitative trend only"
    )
    ax.legend(loc="lower right", framealpha=0.9)

    # Annotate exact values
    for e, te, tl in zip(Ep_uJ, Te_pk, Tl_pk):
        ax.annotate(f"{te:.0f}", (e, te), textcoords="offset points",
                    xytext=(6, 4), fontsize=7, color="#1f77b4")
        ax.annotate(f"{tl:.0f}", (e, tl), textcoords="offset points",
                    xytext=(6, -10), fontsize=7, color="#d62728")

    fig.tight_layout()
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
    cfg = load_config(args.toml)

    print("=" * 78)
    print(f"V1.6 -- 1D TTM prototype  (PROTOTYPE PARAMETERS)")
    print(f"  tau = {cfg.tau_s*1e12:.1f} ps   w0 = {cfg.w0_m*1e6:.1f} um"
          f"   A = {cfg.A}   lambda = {cfg.wavelength_nm:.0f} nm")
    print(f"  delta_abs = {cfg.delta_abs_m*1e9:.1f} nm (placeholder)")
    print(f"  gamma_e = {cfg.gamma_e} J/(m^3 K^2),  ke = {cfg.ke} W/(m K)")
    print(f"  G = {cfg.G:.1e} W/(m^3 K)")
    print(f"  z_max = {cfg.z_max_m*1e9:.0f} nm,  dz = {cfg.dz_m*1e9:.1f} nm")
    print(f"  t in [{cfg.t_start_s*1e12:.0f}, {cfg.t_end_s*1e12:.0f}] ps   "
          f"dt_target = {cfg.dt_s*1e15:.1f} fs")
    print(f"  Ep cases (uJ): {cfg.Ep_uJ_list}")
    print("=" * 78)

    per_case: dict[float, dict] = {}
    z_m_ref: np.ndarray | None = None
    for Ep_uJ, Ep_J in zip(cfg.Ep_uJ_list, cfg.Ep_J_list):
        print(f"\n[run] Ep = {Ep_uJ} uJ ...")
        res = run_ttm_1d(cfg, Ep_J)
        per_case[float(Ep_uJ)] = res
        if z_m_ref is None:
            z_m_ref = res["z_m"]
        Te_pk = float(res["Te_max_K"].max())
        Tl_pk = float(res["Tl_max_K"].max())
        print(f"       -> Te_peak = {Te_pk:.1f} K   Tl_peak = {Tl_pk:.1f} K   "
              f"steps = {res['n_steps']}   dt = {res['dt_actual_s']*1e15:.2f} fs")
        for w in res["warnings"]:
            print(f"       [WARN] {w}", file=sys.stderr)

    # ---- CSV ----
    write_surface_csv(per_case, args.out_root / "ttm_surface_temperature.csv")
    write_depth_csv(per_case, z_m_ref, args.out_root / "ttm_depth_profiles.csv")
    print(f"\n[OK] wrote {args.out_root / 'ttm_surface_temperature.csv'}")
    print(f"[OK] wrote {args.out_root / 'ttm_depth_profiles.csv'}")

    # ---- Figures ----
    fig_dir = args.out_root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    plot_surface_history(per_case, cfg,
                         fig_dir / "v16_ttm_surface_temperature.png")
    plot_depth_profiles(per_case, z_m_ref, cfg,
                        fig_dir / "v16_ttm_depth_profiles.png")
    plot_threshold_map(per_case, cfg,
                       fig_dir / "v16_ttm_threshold_map.png")

    # ---- Headline ----
    print()
    print("TTM headline (PROTOTYPE)")
    print(f"  {'Ep (uJ)':>8}  {'Te_peak (K)':>14}  {'Tl_peak (K)':>14}  "
          f"{'melt?':>5}  {'vapor?':>6}")
    for Ep_uJ, res in per_case.items():
        Te_pk = float(res["Te_max_K"].max())
        Tl_pk = float(res["Tl_max_K"].max())
        melt = "yes" if Tl_pk >= cfg.Tm else "no"
        vapor = "yes" if Tl_pk >= cfg.Tv else "no"
        print(f"  {Ep_uJ:>8g}  {Te_pk:>14.1f}  {Tl_pk:>14.1f}  "
              f"{melt:>5}  {vapor:>6}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
