"""
plot_v62_ttm_md_physical_snapshots.py
=====================================

Create physical model snapshot panels from V6B LAMMPS TTM-MD trajectories.

This is not a workflow diagram.  It plots real MD atom positions from the
LAMMPS trajectory as a center-section physical snapshot.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


KB_EV_K = 8.617333262145e-5
SI_MELT_K = 1687.0
SI_VAP_K = 3538.0


def iter_lammpstrj(path: Path):
    with path.open("r", encoding="utf-8", errors="replace") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(f.readline().strip())
            assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
            natoms = int(f.readline().strip())
            header = f.readline().strip()
            if not header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(header)
            bounds = []
            for _ in range(3):
                lo, hi, *_ = f.readline().split()
                bounds.append((float(lo), float(hi)))
            atom_header = f.readline().strip()
            cols = atom_header.split()[2:]
            data = np.loadtxt([f.readline() for _ in range(natoms)])
            if data.ndim == 1:
                data = data.reshape(1, -1)
            yield step, bounds, cols, data


def parse_times(log_path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    header = None
    with log_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "Step" and parts[1] == "Time":
                header = parts
                continue
            if header and len(parts) >= len(header):
                try:
                    vals = [float(v) for v in parts[: len(header)]]
                except ValueError:
                    header = None
                    continue
                row = dict(zip(header, vals))
                out[int(row["Step"])] = row["Time"]
    return out


def nearest_frames(dump_path: Path, target_steps: list[int]):
    targets = set(target_steps)
    frames = {}
    first = None
    for step, bounds, cols, data in iter_lammpstrj(dump_path):
        if first is None:
            first = (step, bounds, cols, data)
        if step in targets:
            frames[step] = (step, bounds, cols, data)
    if first and first[0] not in frames:
        frames[first[0]] = first
    return frames


def load_summary(summary_path: Path) -> dict[int, dict[str, str]]:
    if not summary_path.exists():
        return {}
    with summary_path.open("r", encoding="utf-8", newline="") as f:
        return {int(r["step"]): r for r in csv.DictReader(f)}


def fmt_time(ps: float) -> str:
    if ps < 1.0:
        return f"{ps:.2f} ps"
    if ps < 100.0:
        return f"{ps:.0f} ps"
    return f"{ps:.0f} ps"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=Path.cwd())
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--case", default="fluence_0p5_100ps")
    ap.add_argument("--section-half-width-A", type=float, default=4.0)
    args = ap.parse_args()

    project_root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else project_root
    case_dir = project_root / "results" / "v61_lammps_ttm_mod_ablation_scan" / args.case
    dump_path = case_dir / f"dump.v61_{args.case}.lammpstrj"
    log_path = case_dir / "log.lammps"
    light_dir = output_root / "lightweight_results" / "v61_lammps_ttm_mod_ablation_scan"
    summary = load_summary(light_dir / f"v61_{args.case}_summary.csv")

    target_times = [5.0, 15.0, 30.0, 60.0, 100.0]
    # V6B dumps every 1 ps after 2 ps equilibration, with dt=0.0001 ps.
    target_steps = [int(round((t + 0.0001) / 0.0001)) for t in target_times]
    # Keep exact known dump steps used in this run.
    target_steps = [50000, 150000, 300000, 600000, 1000000]

    times = parse_times(log_path)
    frames = nearest_frames(dump_path, target_steps)
    if not frames:
        raise FileNotFoundError(dump_path)

    first_step = min(frames)
    _s0, _b0, cols0, data0 = frames[first_step]
    idx0 = {name: i for i, name in enumerate(cols0)}
    initial_surface_A = float(np.nanmax(data0[:, idx0["z"]]))

    fig, axes = plt.subplots(1, len(target_steps), figsize=(18, 5.5), sharey=True, constrained_layout=True)
    fig.suptitle(f"TTM-MD physical snapshots for silicon ablation ({args.case})", fontsize=15)

    for ax, step in zip(axes, target_steps):
        frame = frames.get(step)
        if frame is None:
            ax.axis("off")
            continue
        _step, bounds, cols, data = frame
        idx = {name: i for i, name in enumerate(cols)}
        x = data[:, idx["x"]]
        y = data[:, idx["y"]]
        z = data[:, idx["z"]]
        vz = data[:, idx["vz"]]
        ke = data[:, idx["c_keatom"]]
        temp = np.full_like(ke, np.nan)
        mask_mobile = ke > 1.0e-12
        temp[mask_mobile] = 2.0 * ke[mask_mobile] / (3.0 * KB_EV_K)

        xmid = 0.5 * (bounds[0][0] + bounds[0][1])
        ymid = 0.5 * (bounds[1][0] + bounds[1][1])
        section = np.abs(y - ymid) <= args.section_half_width_A
        signed_x_nm = (x[section] - xmid) * 0.1
        z_nm = (z[section] - initial_surface_A) * 0.1
        t_plot = np.clip(temp[section], 300.0, 8000.0)

        sc = ax.scatter(
            signed_x_nm,
            z_nm,
            c=t_plot,
            s=8,
            cmap="turbo",
            vmin=300.0,
            vmax=8000.0,
            linewidths=0,
            alpha=0.88,
        )
        ejecta = section & (z > initial_surface_A + 10.0) & (vz > 0.0)
        if np.any(ejecta):
            ax.scatter(
                (x[ejecta] - xmid) * 0.1,
                (z[ejecta] - initial_surface_A) * 0.1,
                s=22,
                facecolors="none",
                edgecolors="black",
                linewidths=0.8,
            )

        ax.axhline(0.0, color="white", linewidth=1.0)
        ax.axhline(0.0, color="black", linewidth=0.55, linestyle="--")
        ax.text(0.02, 0.97, fmt_time(times.get(step, math.nan)), transform=ax.transAxes, va="top", ha="left", fontsize=10)

        row = summary.get(step, {})
        ejecta_count = row.get("ejecta_atoms", "0")
        lift = row.get("surface_lift_nm", "")
        ax.text(
            0.02,
            0.04,
            f"ejecta={ejecta_count}\nlift={lift} nm",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.65, "edgecolor": "none", "pad": 2},
        )
        ax.set_xlim(-2.4, 2.4)
        ax.set_ylim(-14.0, 6.0)
        ax.set_xlabel("center section x / nm")
        ax.set_title(f"step {step}", fontsize=9)
        ax.grid(False)
    axes[0].set_ylabel("height from initial surface / nm")

    cbar = fig.colorbar(sc, ax=axes, shrink=0.85, pad=0.015)
    cbar.set_label("atom temperature proxy / K")

    out_dir = output_root / "figures" / "v62_ttm_md_physical_snapshots" / args.case
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"v62_{args.case}_physical_snapshots.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"[OK] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
