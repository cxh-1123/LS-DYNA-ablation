"""
estimate_v16_30ps_scales.py
===========================

V1.6 step 1 -- analytic scale + threshold estimation for the 30 ps branch.

Reads:
  - config/v16_30ps.toml

Writes:
  - results/v16_30ps/v16_30ps_scale_summary.csv
  - results/v16_30ps/v16_30ps_threshold_summary.csv
  - results/v16_30ps/figures/v16_diffusion_length_vs_tau.png
  - results/v16_30ps/figures/v16_threshold_energy_vs_tau.png
  - results/v16_30ps/figures/v16_mesh_resolution_warning.png

Physics:

* Thermal diffusion length:
      L_diff(tau) = sqrt(alpha * tau),  alpha = k / (rho * cp).

* 1D half-infinite body, surface absorption, Gaussian pulse.  The peak
  surface temperature rise for an absorbed surface fluence F_abs is
  (rectangular-pulse approximation, valid within ~1.2x for Gaussian):

      dT_peak = 2 * F_abs / sqrt(pi * tau * rho * cp * k)

  with F_abs = A * F_peak_spatial,  F_peak_spatial = 2 * Ep / (pi * w0^2).
  Inverting for Ep at a target dT:

      Ep(dT) = dT * pi * w0^2 * sqrt(pi * tau * rho * cp * k) / (4 * A)

* "TTM needed" heuristic (silicon-ish):
      tau <  100 ps         -> required
      100 ps <= tau < 1 ns  -> marginal
      tau >= 1 ns           -> not required

Constraints:
  * Pure Python pre-flight, no LS-DYNA, no PyDYNA.
  * Does not modify V1, V1.5, V2, or V3A.
  * Uses only numpy, matplotlib, tomllib.

Run from project root:
    .\\.venv\\Scripts\\python.exe scripts\\estimate_v16_30ps_scales.py
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
from matplotlib.colors import LogNorm


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
    "legend.fontsize": 9,
    "savefig.dpi": 300,
})


# =============================================================================
# Physics helpers (SI internally)
# =============================================================================
def thermal_diffusivity(k: float, rho: float, cp: float) -> float:
    return k / (rho * cp)


def diffusion_length(alpha: float, tau_s: float) -> float:
    return float(np.sqrt(alpha * tau_s))


def Ep_for_dT(dT: float, tau_s: float, w0_m: float, A: float,
              rho: float, cp: float, k: float) -> float:
    """
    Ep in joules required to produce surface peak temperature rise dT on a
    1D half-infinite body with surface absorption and a Gaussian spot.
    """
    return (
        dT * np.pi * w0_m ** 2
        * np.sqrt(np.pi * tau_s * rho * cp * k)
        / (4.0 * A)
    )


def ttm_needed_label(tau_s: float) -> str:
    if tau_s < 100.0e-12:
        return "required (tau < 100 ps)"
    if tau_s < 1.0e-9:
        return "marginal (100 ps <= tau < 1 ns)"
    return "not required (tau >= 1 ns)"


# =============================================================================
# Output writers
# =============================================================================
def write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


# =============================================================================
# Plotting
# =============================================================================
def plot_diffusion_length_vs_tau(
    tau_arr_s: np.ndarray, L_arr_m: np.ndarray,
    dz_rec_arr_m: np.ndarray, v1_dz_m: float,
    discrete_tau_ps: list[float], out_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    tau_ps = tau_arr_s * 1e12
    ax.loglog(tau_ps, L_arr_m * 1e9, color="#1f77b4", linewidth=2.0,
              label="diffusion length  L_diff = sqrt(alpha * tau)")
    ax.loglog(tau_ps, dz_rec_arr_m * 1e9, color="#ff7f0e",
              linewidth=2.0, linestyle="--",
              label="recommended surface dz = L_diff / 5")
    ax.axhline(v1_dz_m * 1e9, color="#d62728", linestyle=":",
               linewidth=1.6, label=f"current V1 dz = {v1_dz_m*1e9:.0f} nm")

    for tau_ps_v in discrete_tau_ps:
        ax.axvline(tau_ps_v, color="0.7", linewidth=0.7, alpha=0.5)
        ax.text(tau_ps_v, ax.get_ylim()[1] * 0.6,
                f" {tau_ps_v:g} ps", rotation=90, va="top", ha="left",
                fontsize=8, color="0.4")

    ax.set_xlabel("pulse width tau (ps)")
    ax.set_ylabel("length (nm)")
    ax.set_title("V1.6 -- thermal diffusion length and recommended dz vs pulse width")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[OK] wrote {out_path}")


def plot_threshold_energy_vs_tau(
    tau_arr_s: np.ndarray,
    Ep_melt_arr_J: np.ndarray, Ep_vap_arr_J: np.ndarray,
    discrete_tau_ps: list[float],
    Ep_scan_uJ: list[float],
    out_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    tau_ps = tau_arr_s * 1e12
    ax.loglog(tau_ps, Ep_melt_arr_J * 1e6,
              color="#1f77b4", linewidth=2.0,
              label="Ep for melt onset (dT = T_m - T0)")
    ax.loglog(tau_ps, Ep_vap_arr_J * 1e6,
              color="#d62728", linewidth=2.0,
              label="Ep for vapor onset (dT = T_v - T0)")

    for tau_ps_v in discrete_tau_ps:
        ax.axvline(tau_ps_v, color="0.7", linewidth=0.7, alpha=0.5)
    for Ep_uJ in Ep_scan_uJ:
        ax.axhline(Ep_uJ, color="0.85", linewidth=0.5, alpha=0.7)
        ax.text(tau_ps[-1] * 0.85, Ep_uJ, f" {Ep_uJ:g} uJ",
                fontsize=7, color="0.4", va="bottom", ha="right")

    ax.set_xlabel("pulse width tau (ps)")
    ax.set_ylabel("required pulse energy Ep (uJ)")
    ax.set_title("V1.6 -- 1D half-infinite analytic threshold Ep vs pulse width\n"
                 "(w0 = 35 um, A = 0.5; surface-absorption approximation)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[OK] wrote {out_path}")


def plot_mesh_resolution_warning(
    tau_arr_s: np.ndarray, dz_rec_arr_m: np.ndarray, v1_dz_m: float,
    out_path: Path,
) -> None:
    """
    A red/green map of under-resolution ratio v1_dz / dz_recommended as
    a function of pulse width.  Anywhere the V1 mesh is more than 1x
    the recommended dz, we are under-resolved -- shaded red.
    """
    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    tau_ps = tau_arr_s * 1e12
    ratio = v1_dz_m / dz_rec_arr_m

    ax.loglog(tau_ps, ratio, color="#9467bd", linewidth=2.0,
              label="V1 dz / recommended dz")
    ax.axhline(1.0, color="#2ca02c", linestyle="--", linewidth=1.5,
               label="adequately resolved (ratio = 1)")
    ax.axhline(2.0, color="#ff7f0e", linestyle="--", linewidth=1.0,
               label="marginal (ratio = 2)")
    ax.fill_between(tau_ps, 1.0, ratio.max() * 1.1, where=(ratio > 1.0),
                    color="#d62728", alpha=0.12,
                    label="V1 mesh under-resolves this pulse width")

    # Annotate ratios at canonical taus
    canonical_tau_ps = [30.0, 100.0, 1000.0, 100000.0]
    for tau_v in canonical_tau_ps:
        tau_s = tau_v * 1e-12
        # interp dz_rec at this tau
        L = np.sqrt((v1_dz_m * 5.0) ** 2)  # placeholder, properly interp below
        idx = int(np.argmin(np.abs(tau_ps - tau_v)))
        r = ratio[idx]
        ax.plot([tau_v], [r], marker="o", color="#9467bd", markersize=7)
        ax.text(tau_v * 1.10, r * 0.9,
                f"x{r:.1f}", fontsize=9, color="#5e3a96")

    ax.set_xlabel("pulse width tau (ps)")
    ax.set_ylabel("under-resolution ratio (V1 dz / recommended dz)")
    ax.set_title("V1.6 -- is the current V1 mesh adequate for this pulse width?\n"
                 "Ratio > 1 means V1 dz is too coarse and a finer local mesh is required")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", framealpha=0.9)

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
    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)

    # Material (single-temperature)
    T0   = float(cfg["material"]["T0_K"])
    Tm   = float(cfg["material"]["T_melt_K"])
    Tv   = float(cfg["material"]["T_vap_K"])
    rho  = float(cfg["material"]["rho_kg_m3"])
    cp   = float(cfg["material"]["cp_J_kgK"])
    k    = float(cfg["material"]["k_W_mK"])
    alpha = thermal_diffusivity(k, rho, cp)

    # Laser
    w0_um = float(cfg["laser"]["spot_radius_um"])
    A     = float(cfg["laser"]["absorption_A"])
    w0_m  = w0_um * 1e-6
    Ep_scan_uJ = list(cfg["laser"]["pulse_energy_uJ"])
    tau_list_ps = list(cfg["laser"]["pulse_widths_ps"])
    cells_per_L = int(cfg["mesh_baseline"]["recommended_cells_per_diffusion_length"])
    v1_dz_um = float(cfg["mesh_baseline"]["current_v1_dz_um"])
    v1_dz_m  = v1_dz_um * 1e-6

    print("=" * 78)
    print("V1.6 -- analytic scale + threshold estimation")
    print(f"  alpha = k/(rho*cp) = {alpha:.4e} m^2/s")
    print(f"  T_m = {Tm} K, T_v = {Tv} K, T0 = {T0} K")
    print(f"  w0 = {w0_um} um,  A = {A}")
    print(f"  V1 baseline dz = {v1_dz_um} um  =  {v1_dz_um*1e3:.0f} nm")
    print(f"  recommended cells per diffusion length: {cells_per_L}")
    print("=" * 78)

    # --- continuous tau grid for plotting ---
    tau_arr_s = np.logspace(np.log10(10e-12), np.log10(1e-3), 200)
    L_arr_m  = np.array([diffusion_length(alpha, t) for t in tau_arr_s])
    dz_rec_m = L_arr_m / cells_per_L

    Ep_melt_arr_J = np.array(
        [Ep_for_dT(Tm - T0, t, w0_m, A, rho, cp, k) for t in tau_arr_s]
    )
    Ep_vap_arr_J = np.array(
        [Ep_for_dT(Tv - T0, t, w0_m, A, rho, cp, k) for t in tau_arr_s]
    )

    # --- discrete tau table ---
    scale_rows: list[dict] = []
    threshold_rows: list[dict] = []
    for tau_ps_v in tau_list_ps:
        tau_s = tau_ps_v * 1e-12
        L_m = diffusion_length(alpha, tau_s)
        dz_rec_m_v = L_m / cells_per_L
        ratio = v1_dz_m / dz_rec_m_v

        scale_rows.append({
            "tau_ps":                f"{tau_ps_v:.3g}",
            "L_diff_um":             f"{L_m*1e6:.4g}",
            "L_diff_nm":             f"{L_m*1e9:.3g}",
            "dz_recommended_nm":     f"{dz_rec_m_v*1e9:.3g}",
            "v1_dz_nm":              f"{v1_dz_m*1e9:.0f}",
            "underres_ratio":        f"{ratio:.2f}",
            "ttm_needed":            ttm_needed_label(tau_s),
            "use_current_v1_mesh":   "no" if ratio > 1.0 else "yes",
        })

        Ep_m = Ep_for_dT(Tm - T0, tau_s, w0_m, A, rho, cp, k)
        Ep_v = Ep_for_dT(Tv - T0, tau_s, w0_m, A, rho, cp, k)
        threshold_rows.append({
            "tau_ps":           f"{tau_ps_v:.3g}",
            "w0_um":            f"{w0_um:.2f}",
            "A":                f"{A:.2f}",
            "Ep_melt_uJ":       f"{Ep_m*1e6:.4g}",
            "Ep_vap_uJ":        f"{Ep_v*1e6:.4g}",
        })

    # --- write CSVs ---
    args.out_root.mkdir(parents=True, exist_ok=True)
    scale_csv     = args.out_root / "v16_30ps_scale_summary.csv"
    threshold_csv = args.out_root / "v16_30ps_threshold_summary.csv"
    write_csv(scale_rows, scale_csv)
    write_csv(threshold_rows, threshold_csv)
    print(f"\n[OK] wrote {scale_csv}")
    print(f"[OK] wrote {threshold_csv}")

    # --- figures ---
    fig_dir = args.out_root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    plot_diffusion_length_vs_tau(
        tau_arr_s, L_arr_m, dz_rec_m, v1_dz_m,
        tau_list_ps,
        fig_dir / "v16_diffusion_length_vs_tau.png",
    )
    plot_threshold_energy_vs_tau(
        tau_arr_s, Ep_melt_arr_J, Ep_vap_arr_J,
        tau_list_ps, Ep_scan_uJ,
        fig_dir / "v16_threshold_energy_vs_tau.png",
    )
    plot_mesh_resolution_warning(
        tau_arr_s, dz_rec_m, v1_dz_m,
        fig_dir / "v16_mesh_resolution_warning.png",
    )

    # --- headline print ---
    print()
    print("Per-tau headline (analytic, 1D half-infinite, w0=35 um, A=0.5)")
    print(f"  {'tau':>8} {'L_diff':>11} {'dz_rec':>10} "
          f"{'V1/rec':>8} {'Ep_melt':>10} {'Ep_vap':>10}  TTM?")
    for s, t in zip(scale_rows, threshold_rows):
        print(
            f"  {s['tau_ps']:>8}  {s['L_diff_nm']:>9} nm "
            f"{s['dz_recommended_nm']:>8} nm "
            f"{s['underres_ratio']:>8} "
            f"{t['Ep_melt_uJ']:>8} uJ "
            f"{t['Ep_vap_uJ']:>8} uJ  "
            f"{s['ttm_needed']}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
