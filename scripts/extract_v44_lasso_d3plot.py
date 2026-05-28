"""
extract_v44_lasso_d3plot.py
===========================

Read real LS-DYNA d3plot fields with lasso-python and export small CSV/PNG
summaries.

Install dependency once:
    python -m pip install lasso-python

Run from project root:
    python scripts\\extract_v44_lasso_d3plot.py --case ep_5p0uj_recoil_5000ms
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import tomllib
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


def import_lasso():
    try:
        from lasso.dyna import ArrayType, D3plot  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on local install
        raise SystemExit(
            "Missing dependency: lasso-python\n"
            "Install it once with:\n"
            "    python -m pip install lasso-python\n"
        ) from exc
    return D3plot, ArrayType


def load_toml(path: Path) -> dict:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def get_array(arrays: dict, array_type, name: str):
    key = getattr(array_type, name, name)
    if key in arrays:
        return arrays[key]
    if name in arrays:
        return arrays[name]
    return None


def selected_case_names(cfg: dict) -> list[str]:
    return [str(row["name"]) for row in cfg.get("selected_cases", [])]


def state_numbers(cfg: dict) -> list[int]:
    return [int(v) for v in cfg["states"]["state_numbers"]]


def shell_centroids(coords: np.ndarray, shell_node_indexes: np.ndarray) -> np.ndarray:
    idx = np.asarray(shell_node_indexes, dtype=np.int64)
    if idx.ndim == 1:
        idx = idx.reshape((-1, 4))
    # Some readers store repeated last nodes for triangular shells; averaging is fine.
    return coords[idx[:, :4], :2].mean(axis=1)


def write_node_displacement_csv(path: Path, node_ids: np.ndarray, coords: np.ndarray,
                                current_or_disp: np.ndarray) -> tuple[float, float, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    # lasso-python exposes the d3plot nodal "displacement" block as current
    # nodal coordinates for these decks.  Convert it to true displacement.
    disp = current_or_disp - coords
    mag = np.linalg.norm(disp, axis=1)
    max_i = int(np.argmax(mag))
    with path.open("w", encoding="utf-8", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["node_id", "x_mm", "y_mm", "z_mm", "ux_mm", "uy_mm", "uz_mm", "u_mag_nm"])
        for nid, xyz, uvw, m in zip(node_ids, coords, disp, mag):
            wr.writerow([
                int(nid),
                f"{xyz[0]:.9e}", f"{xyz[1]:.9e}", f"{xyz[2]:.9e}",
                f"{uvw[0]:.9e}", f"{uvw[1]:.9e}", f"{uvw[2]:.9e}",
                f"{m * 1.0e6:.9e}",
            ])
    return float(mag.min()), float(mag.max()), int(node_ids[max_i])


def write_shell_pressure_csv(path: Path, shell_ids: np.ndarray, centroids: np.ndarray,
                             pressure: np.ndarray) -> tuple[float, float, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if pressure.ndim == 2:
        pressure_plot = pressure.mean(axis=1)
    else:
        pressure_plot = pressure
    max_i = int(np.argmax(pressure_plot))
    with path.open("w", encoding="utf-8", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["element_id", "x_mm", "y_mm", "pressure_GPa"])
        for eid, xy, p in zip(shell_ids, centroids, pressure_plot):
            wr.writerow([int(eid), f"{xy[0]:.9e}", f"{xy[1]:.9e}", f"{p:.9e}"])
    return float(pressure_plot.min()), float(pressure_plot.max()), int(shell_ids[max_i])


def scatter(path: Path, x: np.ndarray, y: np.ndarray, c: np.ndarray,
            title: str, label: str, cmap: str = "turbo") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.5, 3.2), dpi=180)
    sc = ax.scatter(x, y, c=c, s=7, cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel("x / mm")
    ax.set_ylabel("y / mm")
    ax.set_aspect("equal", adjustable="box")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(label)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def pressure_from_shell_stress(shell_stress: np.ndarray) -> np.ndarray:
    # LS-DYNA stress order: xx, yy, zz, xy, yz, xz.
    # Pressure is the negative mean normal stress.
    return -(shell_stress[..., 0] + shell_stress[..., 1] + shell_stress[..., 2]) / 3.0


def extract_case(root: Path, output_root: Path, cfg: dict, case: str) -> list[dict[str, str]]:
    D3plot, ArrayType = import_lasso()
    d3plot_path = root / cfg["source"]["result_root"] / case / "d3plot"
    if not d3plot_path.is_file():
        raise SystemExit(f"d3plot not found: {d3plot_path}")

    states_1based = state_numbers(cfg)
    state_filter = {s - 1 for s in states_1based}
    d3 = D3plot(
        str(d3plot_path),
        state_filter=state_filter,
        state_array_filter=[
            "node_displacement",
            "node_velocity",
            "element_shell_stress",
        ],
    )

    arrays = d3.arrays
    node_ids = get_array(arrays, ArrayType, "node_ids")
    coords = get_array(arrays, ArrayType, "node_coordinates")
    disp = get_array(arrays, ArrayType, "node_displacement")
    vel = get_array(arrays, ArrayType, "node_velocity")
    shell_ids = get_array(arrays, ArrayType, "element_shell_ids")
    shell_idx = get_array(arrays, ArrayType, "element_shell_node_indexes")
    shell_stress = get_array(arrays, ArrayType, "element_shell_stress")
    times = get_array(arrays, ArrayType, "global_timesteps")

    if node_ids is None or coords is None or disp is None:
        raise SystemExit("d3plot did not contain node displacement arrays")

    field_root = output_root / cfg["output"]["field_root"] / case
    fig_root = output_root / cfg["output"]["figure_root"] / case
    rows: list[dict[str, str]] = []

    centroids = None
    if shell_ids is not None and shell_idx is not None:
        centroids = shell_centroids(coords, shell_idx)

    loaded_count = disp.shape[0]
    for local_i, state in enumerate(states_1based[:loaded_count]):
        time_ms = ""
        if times is not None and local_i < len(times):
            time_ms = f"{float(times[local_i]):.9e}"

        state_dir = field_root / f"state_{state:04d}"
        u_csv = state_dir / "node_displacement.csv"
        umin, umax, umax_id = write_node_displacement_csv(u_csv, node_ids, coords, disp[local_i])
        mag_nm = np.linalg.norm(disp[local_i], axis=1) * 1.0e6
        u_fig = fig_root / f"state_{state:04d}_result_displacement_nm.png"
        scatter(u_fig, coords[:, 0], coords[:, 1], mag_nm,
                f"{case} state {state}: real displacement", "resultant displacement / nm")
        rows.append({
            "case": case,
            "state": str(state),
            "time_ms": time_ms,
            "field": "resultant_displacement",
            "min_value": f"{umin:.9e}",
            "max_value": f"{umax:.9e}",
            "unit": "mm",
            "max_id": str(umax_id),
            "csv": str(u_csv),
            "figure": str(u_fig),
        })

        if vel is not None and local_i < vel.shape[0]:
            vmag = np.linalg.norm(vel[local_i], axis=1)
            rows.append({
                "case": case,
                "state": str(state),
                "time_ms": time_ms,
                "field": "resultant_velocity",
                "min_value": f"{float(vmag.min()):.9e}",
                "max_value": f"{float(vmag.max()):.9e}",
                "unit": "mm/ms",
                "max_id": str(int(node_ids[int(np.argmax(vmag))])),
                "csv": "",
                "figure": "",
            })

        if shell_stress is not None and shell_ids is not None and centroids is not None:
            pressure = pressure_from_shell_stress(shell_stress[local_i])
            p_csv = state_dir / "shell_pressure.csv"
            pmin, pmax, pmax_id = write_shell_pressure_csv(p_csv, shell_ids, centroids, pressure)
            p_plot = pressure.mean(axis=1) if pressure.ndim == 2 else pressure
            p_fig = fig_root / f"state_{state:04d}_pressure_GPa.png"
            scatter(p_fig, centroids[:, 0], centroids[:, 1], p_plot,
                    f"{case} state {state}: pressure from real shell stress", "pressure / GPa")
            rows.append({
                "case": case,
                "state": str(state),
                "time_ms": time_ms,
                "field": "pressure_from_shell_stress",
                "min_value": f"{pmin:.9e}",
                "max_value": f"{pmax:.9e}",
                "unit": "GPa",
                "max_id": str(pmax_id),
                "csv": str(p_csv),
                "figure": str(p_fig),
            })
    return rows


def main() -> int:
    default_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=default_root)
    ap.add_argument("--output-root", type=Path, default=None)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--case", default="all", help="case name or all")
    args = ap.parse_args()

    root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else root
    cfg = load_toml(args.config or root / "config" / "v44_lsprepost_real_field_extract.toml")
    requested = selected_case_names(cfg) if args.case == "all" else [args.case]

    rows: list[dict[str, str]] = []
    for case in requested:
        print(f"[READ] {case}")
        rows.extend(extract_case(root, output_root, cfg, case))

    out_dir = output_root / cfg["output"]["lightweight_root"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "v44_real_field_summary.csv"
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = ["case", "state", "time_ms", "field", "min_value", "max_value", "unit", "max_id", "csv", "figure"]
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        wr.writerows(rows)

    (out_dir / "README.md").write_text(
        "# V44 real d3plot field summary\n\n"
        "Generated by `scripts/extract_v44_lasso_d3plot.py` from real d3plot arrays.\n\n"
        f"- rows: {len(rows)}\n"
        "- displacement comes from `node_displacement`\n"
        "- pressure is computed as `-(sxx + syy + szz) / 3` from `element_shell_stress`\n",
        encoding="utf-8",
    )
    print(f"[OK] {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
