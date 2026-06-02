"""
plot_v61_ttm_mod_evolution.py
=============================

Plot V6B LAMMPS ttm/mod ablation scout summary.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def val(row: dict[str, str], key: str) -> float:
    txt = row.get(key, "")
    return float(txt) if txt else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=Path.cwd())
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="fluence_0p5_100ps")
    args = ap.parse_args()
    project_root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else project_root
    light_dir = output_root / "lightweight_results" / "v61_lammps_ttm_mod_ablation_scan"
    fig_dir = output_root / "figures" / "v61_lammps_ttm_mod_ablation_scan" / args.case
    fig_dir.mkdir(parents=True, exist_ok=True)

    data = rows(light_dir / f"v61_{args.case}_summary.csv")
    t = [val(r, "time_ps") for r in data]
    thermo_t = [val(r, "thermo_temp_K") for r in data]
    atom_tmax = [val(r, "max_atom_temperature_K") for r in data]
    etmax = [val(r, "max_electron_temperature_K") for r in data]
    etmean = [val(r, "mean_electron_temperature_K") for r in data]
    lift = [val(r, "surface_lift_nm") for r in data]
    ejecta = [val(r, "ejecta_atoms") for r in data]
    melt = [val(r, "melt_like_atom_count") for r in data]
    vap = [val(r, "vap_like_atom_count") for r in data]
    spmax = [val(r, "max_stress_proxy_eV") for r in data]
    spmin = [val(r, "min_stress_proxy_eV") for r in data]

    fig, axs = plt.subplots(2, 3, figsize=(16, 8), constrained_layout=True)
    fig.suptitle(f"V6B ttm/mod Si ablation scout: {args.case}", fontsize=16)

    ax = axs[0, 0]
    ax.plot(t, etmax, label="max electron T", color="#d62728")
    ax.plot(t, etmean, label="mean electron T", color="#ff7f0e")
    ax.plot(t, thermo_t, label="lattice T (LAMMPS thermo)", color="#1f77b4")
    ax.axhline(1687, color="k", linestyle="--", linewidth=1, label="Si melt")
    ax.axhline(3538, color="0.4", linestyle=":", linewidth=1, label="Si vapor")
    ax.set_xlabel("time / ps")
    ax.set_ylabel("temperature / K")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axs[0, 1]
    ax.plot(t, atom_tmax, color="#17becf", label="max atom T proxy")
    ax.axhline(1687, color="k", linestyle="--", linewidth=1)
    ax.axhline(3538, color="0.4", linestyle=":", linewidth=1)
    ax.set_xlabel("time / ps")
    ax.set_ylabel("temperature / K")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axs[0, 2]
    ax.plot(t, melt, label="T > melt", color="#ff7f0e")
    ax.plot(t, vap, label="T > vapor", color="#d62728")
    ax.set_xlabel("time / ps")
    ax.set_ylabel("atom count")
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

    ax = axs[1, 2]
    ax.plot(t, spmax, label="max stress proxy", color="#d62728")
    ax.plot(t, spmin, label="min stress proxy", color="#1f77b4")
    ax.axhline(0.0, color="k", linewidth=0.8)
    ax.set_xlabel("time / ps")
    ax.set_ylabel("stress proxy / eV")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    out = fig_dir / f"v61_{args.case}_overview.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
