"""
plot_v17_30ps_temperature.py
============================

Post-process V1.7 local 30 ps LS-DYNA temperature fields from tprint.

Outputs (300 dpi PNG, English labels):
  results/v17_30ps_local/figures/v17_Tmax_vs_time.png
  results/v17_30ps_local/figures/v17_temperature_snapshots.png
  results/v17_30ps_local/figures/v17_full_section_30ps.png
  results/v17_30ps_local/figures/v17_full_section_peak.png
  results/v17_30ps_local/figures/v17_zoomed_center_hot_zone.png
  results/v17_30ps_local/figures/v17_zoomed_center_hot_zone_peak.png
  results/v17_30ps_local/figures/v17_surface_radial_profiles.png
  results/v17_30ps_local/figures/v17_depth_profiles.png

Requires LS-DYNA results in results/v17_30ps_local/<case>/tprint.
If no results exist, prints a clear message and exits non-zero.

Run:
    .\\.venv\\Scripts\\python.exe scripts\\plot_v17_30ps_temperature.py
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
from matplotlib.colors import PowerNorm

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from plot_v1_temperature import parse_tprint  # noqa: E402
from check_v17_outputs import load_mesh_nodes, tprint_to_grid  # noqa: E402


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

CASE_STYLE = {
    "ep_0p5uj": {"color": "#1f77b4", "marker": "o"},
    "ep_1p0uj": {"color": "#2ca02c", "marker": "s"},
    "ep_2p0uj": {"color": "#ff7f0e", "marker": "^"},
    "ep_5p0uj": {"color": "#d62728", "marker": "D"},
}


def load_mesh_arrays(mesh_csv: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    r_um, z_um = [], []
    node_map: dict[int, tuple[int, int]] = {}
    with mesh_csv.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            i, j = int(row["i"]), int(row["j"])
            node_map[int(row["node_id"])] = (i, j)
            if len(r_um) <= i:
                r_um.append(float(row["r_um"]))
            if len(z_um) <= j:
                z_um.append(float(row["z_um"]))
    return np.asarray(r_um), np.asarray(z_um), node_map


def mirror_field(r_um: np.ndarray, T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r_full = np.concatenate([-r_um[:0:-1], r_um])
    T_full = np.concatenate([T[:, :0:-1], T], axis=1)
    return r_full, T_full


def load_case_data(case_dir: Path, node_map, NR, NZ):
    tprint = case_dir / "tprint"
    if not tprint.is_file():
        return None
    times_ms, T_all = parse_tprint(tprint)
    times_ns = np.asarray(times_ms) * 1e6
    grids = [tprint_to_grid(T_all[k], node_map, NR, NZ) for k in range(T_all.shape[0])]
    Tmax = np.array([float(g.max()) for g in grids])
    return {"times_ns": times_ns, "grids": grids, "Tmax": Tmax}


def print_case_time_diagnostics(case_name: str, times_ns: np.ndarray) -> None:
    """Print tprint time coverage for one case."""
    times_ns = np.asarray(times_ns, dtype=float)
    n = len(times_ns)
    head = times_ns[: min(10, n)]
    print(
        f"  [{case_name}] snapshots={n}, "
        f"min_time_ns={float(times_ns[0]):.6g}, "
        f"max_time_ns={float(times_ns[-1]):.6g}"
    )
    print(f"    first {len(head)} times (ns): {[float(t) for t in head]}")


def pick_nearest_time_index(
    times_ns: np.ndarray,
    target_ns: float,
    *,
    case_name: str = "",
) -> tuple[int, float, float]:
    """
    Return (index, actual_time_ns, abs_error_ns) for the snapshot nearest target_ns.

    - target_ns and times_ns are both in nanoseconds.
    - If target is outside [min, max], clamp to the boundary frame and warn.
    - For target > 0, skip the t=0 initial frame when later frames exist so
      sparse tprint (first dump after t=0 often >> 0.03 ns) does not stick at 300 K.
    - target == 0 explicitly selects the initial frame.
    """
    times_ns = np.asarray(times_ns, dtype=float)
    if times_ns.size == 0:
        raise ValueError("times_ns is empty")

    t_min = float(times_ns[0])
    t_max = float(times_ns[-1])
    label = f" ({case_name})" if case_name else ""

    if target_ns > t_max:
        print(
            f"  [WARN] target time {target_ns:g} ns outside available range "
            f"[{t_min:g}, {t_max:g}] ns{label}; using last frame",
            file=sys.stderr,
        )
        idx = int(times_ns.size - 1)
        actual = float(times_ns[idx])
        return idx, actual, abs(actual - target_ns)

    if target_ns < t_min:
        print(
            f"  [WARN] target time {target_ns:g} ns outside available range "
            f"[{t_min:g}, {t_max:g}] ns{label}; using first frame",
            file=sys.stderr,
        )
        actual = t_min
        return 0, actual, abs(actual - target_ns)

    if target_ns == 0.0:
        return 0, t_min, abs(t_min - target_ns)

    pool = np.arange(times_ns.size, dtype=int)
    if t_min == 0.0 and times_ns.size > 1:
        pool = pool[1:]

    diffs = np.abs(times_ns[pool] - target_ns)
    idx = int(pool[int(np.argmin(diffs))])
    actual = float(times_ns[idx])
    return idx, actual, abs(actual - target_ns)


def pick_peak_time_index(times_ns: np.ndarray, Tmax: np.ndarray) -> tuple[int, float]:
    """Return (index, actual_time_ns) at global Tmax peak."""
    idx = int(np.argmax(Tmax))
    return idx, float(times_ns[idx])


def format_target_actual_title(
    target_ns: float,
    actual_ns: float,
    *,
    extra: str = "",
) -> str:
    parts = [f"target={target_ns:g} ns", f"actual={actual_ns:g} ns"]
    if extra:
        parts.append(extra)
    return ", ".join(parts)


def format_time_legend_label(actual_ns: float) -> str:
    """Human-readable time for legend (ps if < 1 ns, else ns)."""
    if actual_ns < 1.0:
        return f"t={actual_ns * 1e3:.2f} ps"
    return f"t={actual_ns:g} ns"


def get_output_times_ns(cfg: dict) -> list[float]:
    """
    Read output snapshot times (ns) from the V1.7 TOML config.

    Supports both layouts:

        [output_times_ns]
        times = [0.0, 0.01, ...]

    and:

        output_times_ns = [0.0, 0.01, ...]
    """
    raw = cfg.get("output_times_ns")
    if raw is None:
        raise KeyError(
            'config key "output_times_ns" is missing.  '
            "Add [output_times_ns] times = [...] or output_times_ns = [...]."
        )

    if isinstance(raw, list):
        times_raw = raw
    elif isinstance(raw, dict):
        if "times" in raw:
            times_raw = raw["times"]
        else:
            raise ValueError(
                f'output_times_ns is a dict but has no "times" key.  '
                f"type={type(raw).__name__}, value={raw!r}, keys={list(raw.keys())}"
            )
    else:
        raise ValueError(
            f"output_times_ns has unsupported type {type(raw).__name__}: {raw!r}.  "
            "Expected a list of floats or [output_times_ns] with times = [...]."
        )

    if not isinstance(times_raw, list):
        raise ValueError(
            f"output_times_ns resolved to non-list type {type(times_raw).__name__}: "
            f"{times_raw!r}"
        )

    out: list[float] = []
    for item in times_raw:
        try:
            out.append(float(item))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"output_times_ns contains non-numeric entry {item!r} "
                f"(type {type(item).__name__}).  Full list: {times_raw!r}"
            ) from exc
    return out


def draw_not_run_panel(ax, case_name: str) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.5, 0.5,
        f"{case_name}\n(not run -- no tprint)",
        ha="center", va="center", fontsize=11, color="0.35",
        transform=ax.transAxes,
    )
    ax.set_title(case_name, fontsize=10, color="0.45")


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--toml", type=Path,
                    default=project_root / "config" / "v17_30ps_local_mesh.toml")
    ap.add_argument("--registry", type=Path,
                    default=project_root / "models" / "v17_30ps_local" / "v17_case_registry.csv")
    ap.add_argument("--mesh-nodes", type=Path,
                    default=project_root / "models" / "v17_30ps_local" / "v17_mesh_nodes.csv")
    ap.add_argument("--out-root", type=Path,
                    default=project_root / "results" / "v17_30ps_local")
    args = ap.parse_args()

    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)
    T_melt = float(cfg["material"]["T_melt_K"])
    T_vap = float(cfg["material"]["T_vap_K"])
    output_times_ns = get_output_times_ns(cfg)
    target_30ps_ns = 0.03

    with args.registry.open("r", encoding="utf-8") as fh:
        registry = list(csv.DictReader(fh))

    scheduled = [
        row for row in registry
        if row["run_lsdyna_final"].strip().lower() in ("true", "1", "yes")
    ]
    scheduled_names = [row["name"] for row in scheduled]

    r_um, z_um, node_map = load_mesh_arrays(args.mesh_nodes)
    NR = len(r_um) - 1
    NZ = len(z_um) - 1

    cases: dict[str, dict] = {}
    missing_cases: list[str] = []
    for row in scheduled:
        name = row["name"]
        d = load_case_data(args.out_root / name, node_map, NR, NZ)
        if d is not None:
            cases[name] = d
        else:
            missing_cases.append(name)
            print(f"  [SKIP] {name}: no tprint -- panel will show 'not run'", file=sys.stderr)

    if not cases:
        print("[ERROR] No V1.7 tprint files found under results/v17_30ps_local/.",
              file=sys.stderr)
        print("        Run scripts\\run_v17_selected.ps1 first.", file=sys.stderr)
        return 3

    if missing_cases:
        print(f"  [INFO] plotting {len(cases)} case(s); "
              f"{len(missing_cases)} missing: {', '.join(missing_cases)}")

    print("\n[tprint time diagnostics]")
    for name, d in cases.items():
        print_case_time_diagnostics(name, d["times_ns"])

    fig_dir = args.out_root / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Tmax vs time ----
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for name, d in cases.items():
        st = CASE_STYLE.get(name, {"color": "0.3", "marker": "o"})
        ax.plot(d["times_ns"], d["Tmax"], color=st["color"], marker=st["marker"],
                markevery=max(1, len(d["times_ns"]) // 15), linewidth=1.6,
                label=f"{name} (Ep from registry)")
    ax.axhline(T_melt, color="#00aaff", linestyle="--", linewidth=1.0, label=f"T_m={T_melt:.0f}K")
    ax.axhline(T_vap, color="#ff5252", linestyle="--", linewidth=1.0, label=f"T_v={T_vap:.0f}K")
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("Tmax (K)")
    ax.set_title("V1.7 -- Tmax vs time (single-temperature LS-DYNA local mesh)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    p1 = fig_dir / "v17_Tmax_vs_time.png"
    fig.savefig(p1)
    plt.close(fig)
    print(f"[OK] {p1}")

    # Pick reference case: highest Tmax peak
    ref_name = max(cases.keys(), key=lambda n: cases[n]["Tmax"].max())
    ref = cases[ref_name]
    k30, t30_actual_ns, err30 = pick_nearest_time_index(
        ref["times_ns"], target_30ps_ns, case_name=ref_name,
    )
    T_ref = ref["grids"][k30]
    k_peak_ref, t_peak_ref_ns = pick_peak_time_index(ref["times_ns"], ref["Tmax"])
    T_ref_peak = ref["grids"][k_peak_ref]
    print(
        f"\n[snapshot selection] ref={ref_name}: "
        f"30ps target={target_30ps_ns:g} ns -> actual={t30_actual_ns:g} ns "
        f"(|err|={err30:g} ns); peak at {t_peak_ref_ns:g} ns, Tmax={ref['Tmax'][k_peak_ref]:.1f} K"
    )

    def plot_full_section_pair(
        T_grid: np.ndarray,
        title_suffix: str,
        out_path: Path,
        norm: PowerNorm,
    ) -> None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].pcolormesh(r_um, z_um, T_grid, cmap="inferno", norm=norm, shading="gouraud")
        axes[0].set_xlim(0, 100)
        axes[0].set_ylim(0, 5)
        axes[0].set_xlabel("r (um)")
        axes[0].set_ylabel("z (um)")
        axes[0].set_title(f"Raw half section ({ref_name})\nlaser at upper-left")
        r_full_l, T_full_l = mirror_field(r_um, T_grid)
        im = axes[1].pcolormesh(r_full_l, z_um, T_full_l, cmap="inferno", norm=norm, shading="gouraud")
        axes[1].set_xlim(-100, 100)
        axes[1].set_ylim(0, 5)
        axes[1].set_xlabel("r (um)")
        axes[1].set_ylabel("z (um)")
        axes[1].set_title("Mirrored full section\ncentre hot zone at top centre")
        axes[1].axvline(0, color="cyan", ls="--", lw=0.8)
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, label="T (K)")
        fig.suptitle(f"V1.7 -- full local section ({title_suffix})", y=1.02)
        fig.tight_layout()
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        print(f"[OK] {out_path}")

    def plot_zoomed_pair(
        T_grid: np.ndarray,
        title_suffix: str,
        out_path: Path,
    ) -> None:
        r_full_l, T_full_l = mirror_field(r_um, T_grid)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        norm_z = PowerNorm(
            gamma=0.40,
            vmin=max(300, T_grid[T_grid > 300].min() if (T_grid > 300).any() else 300),
            vmax=T_grid.max(),
        )
        axes[0].pcolormesh(r_um, z_um, T_grid, cmap="inferno", norm=norm_z, shading="gouraud")
        axes[0].set_xlim(0, 50)
        axes[0].set_ylim(4.5, 5.0)
        axes[0].set_title("Zoom raw half: r=0--50 um, z=4.5--5 um")
        axes[1].pcolormesh(r_full_l, z_um, T_full_l, cmap="inferno", norm=norm_z, shading="gouraud")
        axes[1].set_xlim(-50, 50)
        axes[1].set_ylim(4.5, 5.0)
        axes[1].axvline(0, color="cyan", ls="--", lw=0.8)
        axes[1].set_title("Zoom mirrored: r=-50--50 um, z=4.5--5 um\n(zoomed color range)")
        for ax in axes:
            ax.set_xlabel("r (um)")
            ax.set_ylabel("z (um)")
        fig.suptitle(f"V1.7 -- zoomed centre hot zone ({ref_name}, {title_suffix})", y=1.02)
        fig.tight_layout()
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        print(f"[OK] {out_path}")

    # ---- 2. snapshots multi-case at ~30 ps (include "not run" placeholders) ----
    n = max(len(scheduled_names), 1)
    ncols = 2
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4.5 * nrows), squeeze=False)
    Tmax_all = max(
        d["grids"][pick_nearest_time_index(d["times_ns"], target_30ps_ns, case_name=n)[0]].max()
        for n, d in cases.items()
    )
    norm = PowerNorm(gamma=0.45, vmin=300.0, vmax=Tmax_all)
    for ax, name in zip(axes.ravel(), scheduled_names):
        if name not in cases:
            draw_not_run_panel(ax, name)
            continue
        d = cases[name]
        k, actual_ns, _ = pick_nearest_time_index(
            d["times_ns"], target_30ps_ns, case_name=name,
        )
        Tg = d["grids"][k]
        r_full, T_full = mirror_field(r_um, Tg)
        ax.pcolormesh(r_full, z_um, T_full, cmap="inferno", norm=norm, shading="gouraud")
        ax.set_xlim(-100, 100)
        ax.set_ylim(0, 5)
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        title = format_target_actual_title(
            target_30ps_ns, actual_ns, extra=f"Tmax={Tg.max():.1f} K",
        )
        ax.set_title(f"{name}\n{title}")
        ax.axvline(0, color="cyan", ls="--", lw=0.8, alpha=0.7)
    for ax in axes.ravel()[len(scheduled_names):]:
        ax.axis("off")
    fig.suptitle(
        f"V1.7 -- temperature snapshots (target={target_30ps_ns:g} ns, mirrored full section)",
        y=1.0,
    )
    fig.tight_layout()
    p2 = fig_dir / "v17_temperature_snapshots.png"
    fig.savefig(p2, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {p2}")

    # ---- 3. full section at 30 ps (ref case) ----
    title_30 = format_target_actual_title(
        target_30ps_ns, t30_actual_ns,
        extra=f"Tmax={T_ref.max():.1f} K (nearest frame)",
    )
    plot_full_section_pair(
        T_ref, title_30, fig_dir / "v17_full_section_30ps.png", norm,
    )

    # ---- 3b. full section at Tmax peak (ref case) ----
    norm_peak = PowerNorm(gamma=0.45, vmin=300.0, vmax=float(T_ref_peak.max()))
    title_peak = format_target_actual_title(
        t_peak_ref_ns, t_peak_ref_ns,
        extra=f"Tmax={T_ref_peak.max():.1f} K (peak frame)",
    )
    plot_full_section_pair(
        T_ref_peak, title_peak, fig_dir / "v17_full_section_peak.png", norm_peak,
    )

    # ---- 4. zoomed center hot zone at 30 ps ----
    title_z30 = format_target_actual_title(target_30ps_ns, t30_actual_ns)
    plot_zoomed_pair(T_ref, title_z30, fig_dir / "v17_zoomed_center_hot_zone.png")

    # ---- 4b. zoomed center hot zone at peak ----
    title_zpk = format_target_actual_title(
        t_peak_ref_ns, t_peak_ref_ns,
        extra=f"Tmax={T_ref_peak.max():.1f} K",
    )
    plot_zoomed_pair(T_ref_peak, title_zpk, fig_dir / "v17_zoomed_center_hot_zone_peak.png")

    # ---- 5. surface radial profiles ----
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for name, d in cases.items():
        st = CASE_STYLE.get(name, {"color": "0.3", "marker": "o"})
        for target_ns in output_times_ns:
            k, actual_ns, _ = pick_nearest_time_index(
                d["times_ns"], target_ns, case_name=name,
            )
            Ttop = d["grids"][k][-1, :]
            ax.plot(
                r_um, Ttop, color=st["color"], linewidth=1.0, alpha=0.85,
                label=f"{name} {format_time_legend_label(actual_ns)}",
            )
    ax.axhline(T_melt, color="#00aaff", ls="--", lw=1.0)
    ax.axhline(T_vap, color="#ff5252", ls="--", lw=1.0)
    ax.set_xlabel("r (um)  -- top surface")
    ax.set_ylabel("T (K)")
    ax.set_xlim(0, 100)
    ax.set_title("V1.7 -- surface radial temperature profiles")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    p5 = fig_dir / "v17_surface_radial_profiles.png"
    fig.savefig(p5)
    plt.close(fig)
    print(f"[OK] {p5}")

    # ---- 6. depth profiles on axis at config output times ----
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for name, d in cases.items():
        st = CASE_STYLE.get(name, {"color": "0.3", "marker": "o"})
        for target_ns in output_times_ns:
            k, actual_ns, _ = pick_nearest_time_index(
                d["times_ns"], target_ns, case_name=name,
            )
            Tax = d["grids"][k][:, 0]
            ax.plot(
                z_um, Tax, color=st["color"], linewidth=1.6,
                label=f"{name} {format_time_legend_label(actual_ns)}",
            )
    ax.axhline(T_melt, color="#00aaff", ls="--", lw=1.0)
    ax.axhline(T_vap, color="#ff5252", ls="--", lw=1.0)
    ax.set_xlabel("z (um)  -- depth from bottom")
    ax.set_ylabel("T (K) at r=0")
    ax.set_title("V1.7 -- on-axis depth temperature profiles (config output times)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    p6 = fig_dir / "v17_depth_profiles.png"
    fig.savefig(p6)
    plt.close(fig)
    print(f"[OK] {p6}")

    print(f"\n[DONE] figures in {fig_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
