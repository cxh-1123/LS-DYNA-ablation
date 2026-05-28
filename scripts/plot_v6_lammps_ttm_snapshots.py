"""
plot_v6_lammps_ttm_snapshots.py
================================

Plot V6 LAMMPS TTM-MD snapshot trajectories (side-view r-z, T proxy from speed/ke).

Run after LAMMPS completes:
    python scripts\\plot_v6_lammps_ttm_snapshots.py --case ep_5p0uj
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

R_LIM = (-60.0, 60.0)
Z_LIM = (-8.0, 80.0)
MAIN_TITLE = "V6 LAMMPS TTM-MD snapshots (ep_5p0uj, early stage)"


def read_lammpstrj(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Minimal lammpstrj reader -> x,y,z in Angstrom."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    natoms = 0
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    mode = None
    for line in lines:
        if line.strip() == "ITEM: NUMBER OF ATOMS":
            mode = "natoms"
            continue
        if line.startswith("ITEM: ATOMS"):
            mode = "atoms"
            cols = line.split()[2:]
            ix = cols.index("x") if "x" in cols else 2
            iy = cols.index("y") if "y" in cols else 3
            iz = cols.index("z") if "z" in cols else 4
            continue
        if mode == "natoms":
            natoms = int(line.strip())
            mode = None
            continue
        if mode == "atoms":
            parts = line.split()
            if len(parts) < 4:
                continue
            xs.append(float(parts[ix]))
            ys.append(float(parts[iy]))
            zs.append(float(parts[iz]))
            if len(xs) >= natoms:
                break
    return np.asarray(xs), np.asarray(ys), np.asarray(zs)


def find_snapshots(res_dir: Path) -> list[tuple[float, Path]]:
    out: list[tuple[float, Path]] = []
    for p in sorted(res_dir.glob("snapshot_*ps.lammpstrj")):
        m = re.search(r"snapshot_([0-9.]+)ps", p.name)
        if m:
            out.append((float(m.group(1)), p))
    if not out and (res_dir / "traj.lammpstrj").is_file():
        out.append((0.0, res_dir / "traj.lammpstrj"))
    return sorted(out, key=lambda t: t[0])


def _view_coords(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, laser_axis: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if laser_axis == "+x":
        r_um = np.sqrt(y ** 2 + z ** 2) * 0.1
        z_um = -x * 0.1
        speed = np.sqrt(x ** 2 + y ** 2 + z ** 2) * 0.01
    else:
        r_um = x * 0.1
        z_um = z * 0.1
        speed = np.sqrt(x ** 2 + z ** 2) * 0.01
    return r_um, z_um, speed


def plot_sequence(
    snaps: list[tuple[float, Path]], out_path: Path, *, laser_axis: str = "+x"
) -> None:
    avail_times = sorted({round(t, 3) for t, _ in snaps})
    # Pad to 7 panels: use available times, then blanks
    target = [0, 5, 10, 30, 50, 100, 200, 500]
    times: list[float] = []
    for t in target:
        if t in avail_times and t not in times:
            times.append(t)
    for t in avail_times:
        if t not in times and len(times) < 7:
            times.append(t)
    while len(times) < 7:
        times.append(float("nan"))
    times = times[:7]
    snap_map = {round(t, 3): p for t, p in snaps}
    fig, axes = plt.subplots(2, 4, figsize=(14.0, 7.2), dpi=300)
    vmax_ke = 1.0
    for _, p in snaps:
        x, y, z = read_lammpstrj(p)
        _, _, speed = _view_coords(x, y, z, laser_axis)
        vmax_ke = max(vmax_ke, float(speed.max()))

    for ax, t in zip(axes.ravel()[:7], times):
        p = snap_map.get(float(t)) if not np.isnan(t) else None
        ax.set_facecolor("#f4f4f4")
        if p is None:
            label = "not run" if np.isnan(t) else f"no {int(t)} ps snap"
            ax.text(0.5, 0.5, label, transform=ax.transAxes, ha="center")
        else:
            x, y, z = read_lammpstrj(p)
            r_um, z_um, speed = _view_coords(x, y, z, laser_axis)
            ax.scatter(r_um, z_um, c=speed, s=1.5, cmap="inferno", vmin=0, vmax=vmax_ke, alpha=0.65, linewidths=0)
            ax.axhline(0, color="#333", lw=0.8)
        ax.set_xlim(R_LIM)
        ax.set_ylim(Z_LIM)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{int(t)} ps" if not np.isnan(t) else "-", fontsize=9)
        ax.tick_params(labelsize=7)
    axes.ravel()[7].axis("off")
    axes.ravel()[7].text(
        0.5, 0.5,
        "V6 TTM-MD atom snapshots\nlaser axis +x -> height\nRequires LAMMPS run",
        ha="center", va="center", fontsize=8,
    )
    for i in [0, 4]:
        axes.ravel()[i].set_ylabel("z / um")
    for ax in axes.ravel()[:7]:
        ax.set_xlabel("r / um")
    fig.suptitle(MAIN_TITLE, fontsize=12)
    fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.10, hspace=0.25, wspace=0.15)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {out_path}")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", default="ep_5p0uj")
    args = ap.parse_args()

    res_dir = root / "results/v6_lammps_ttm_md_pilot" / args.case
    case_dir = root / "models/v6_lammps_ttm_md_pilot" / args.case
    laser_axis = "+x"
    manifest = case_dir / "v6_run_manifest.json"
    if manifest.is_file():
        meta = json.loads(manifest.read_text(encoding="utf-8"))
        laser_axis = str(meta.get("laser_axis", laser_axis))

    if not res_dir.is_dir():
        print(f"[WARN] no results yet: {res_dir}\n       run scripts/run_v6_selected.ps1 first", file=sys.stderr)
        return 2

    snaps = find_snapshots(res_dir)
    if not snaps:
        snaps = find_snapshots(case_dir)
    if not snaps:
        print(f"[WARN] no snapshot_*.lammpstrj in {res_dir}", file=sys.stderr)
        return 2

    out = root / "figures/v6_lammps_ttm_md_pilot" / args.case / "v6_ttm_md_snapshot_sequence.png"
    plot_sequence(snaps, out, laser_axis=laser_axis)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
