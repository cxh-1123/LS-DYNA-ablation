"""
plot_v60_ttm_md_evolution.py
============================

Plot V6A lightweight LAMMPS TTM-MD summary in the existing paper-style layout.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str) -> float:
    txt = row.get(key, "")
    return float(txt) if txt else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=Path.cwd())
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="ttm_initial_pulse_5uj_equiv")
    args = ap.parse_args()

    project_root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else project_root
    light_dir = output_root / "lightweight_results" / "v60_lammps_ttm_md_pilot"
    fig_dir = output_root / "figures" / "v60_lammps_ttm_md_pilot" / args.case
    fig_dir.mkdir(parents=True, exist_ok=True)

    rows = read_rows(light_dir / "v60_ttm_md_summary.csv")
    rows = [r for r in rows if r["case"] == args.case]

    t = [f(r, "time_ps") for r in rows]
    atom_t = [f(r, "mean_atom_temperature_K") for r in rows]
    atom_t_max = [f(r, "max_atom_temperature_K") for r in rows]
    elec_t = [f(r, "mean_electron_temperature_K") for r in rows]
    elec_t_max = [f(r, "max_electron_temperature_K") for r in rows]
    pmax = [f(r, "max_stress_proxy_eV") for r in rows]
    pmin = [f(r, "min_stress_proxy_eV") for r in rows]
    lift = [f(r, "surface_lift_nm") for r in rows]
    ejecta = [f(r, "ejecta_atoms") for r in rows]

    fig, axs = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
    fig.suptitle("V6A TTM-MD silicon pilot: electron-lattice coupling check", fontsize=16)

    ax = axs[0, 0]
    ax.plot(t, elec_t, label="mean electron T", color="#d62728")
    ax.plot(t, elec_t_max, label="max electron T", color="#ff7f0e", linestyle="--")
    ax.plot(t, atom_t, label="mean lattice T", color="#1f77b4")
    ax.plot(t, atom_t_max, label="max lattice T", color="#17becf", linestyle="--")
    ax.set_xlabel("time / ps")
    ax.set_ylabel("temperature / K")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axs[0, 1]
    ax.plot(t, pmax, label="max stress proxy", color="#d62728")
    ax.plot(t, pmin, label="min stress proxy", color="#1f77b4")
    ax.axhline(0.0, color="k", linewidth=0.8)
    ax.set_xlabel("time / ps")
    ax.set_ylabel("stress proxy / eV")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axs[1, 0]
    ax.plot(t, lift, color="#2ca02c")
    ax.set_xlabel("time / ps")
    ax.set_ylabel("surface lift / nm")
    ax.grid(True, alpha=0.3)

    ax = axs[1, 1]
    ax.step(t, ejecta, where="post", color="#9467bd")
    ax.set_xlabel("time / ps")
    ax.set_ylabel("ejecta atoms")
    ax.grid(True, alpha=0.3)

    out = fig_dir / "v60_ttm_md_overview.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
