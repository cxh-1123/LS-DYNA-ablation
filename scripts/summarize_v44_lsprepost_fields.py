"""
summarize_v44_lsprepost_fields.py
=================================

Summarize text field files exported by LS-PrePost V44 command files.

Run after LS-PrePost export:
    python scripts\\summarize_v44_lsprepost_fields.py
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
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def parse_keyword_mesh(path: Path) -> tuple[dict[int, tuple[float, float, float]], dict[int, tuple[int, ...]]]:
    nodes: dict[int, tuple[float, float, float]] = {}
    shells: dict[int, tuple[int, ...]] = {}
    section = ""
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line or line.startswith("$"):
            continue
        if line.startswith("*"):
            up = line.upper()
            if up.startswith("*NODE"):
                section = "node"
            elif up.startswith("*ELEMENT_SHELL"):
                section = "shell"
            else:
                section = ""
            continue
        parts = line.replace(",", " ").split()
        if section == "node" and len(parts) >= 4:
            nodes[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))
        elif section == "shell" and len(parts) >= 6:
            eid = int(parts[0])
            nids = tuple(int(v) for v in parts[2:6])
            shells[eid] = nids
    return nodes, shells


def parse_nodal_displacement(path: Path) -> dict[int, tuple[float, float, float]]:
    out: dict[int, tuple[float, float, float]] = {}
    if not path.is_file():
        return out
    active = False
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if line.upper().startswith("$NODAL_DISPLACEMENT"):
            active = True
            continue
        if not active or not line or line.startswith("$") or line.startswith("*"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) >= 4:
            out[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))
    return out


def parse_scalar_results(path: Path) -> tuple[str, dict[int, float]]:
    name = path.stem
    out: dict[int, float] = {}
    if not path.is_file():
        return name, out
    active = False
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if line.startswith("$RESULT OF"):
            name = line.replace("$RESULT OF", "").strip() or name
            continue
        if line.upper().startswith(("$NODAL_RESULTS", "$SHELL_ELEMENT_RESULTS")):
            active = True
            continue
        if not active or not line or line.startswith("$") or line.startswith("*"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) >= 2:
            out[int(parts[0])] = float(parts[1])
    return name, out


def case_names(cfg: dict) -> list[str]:
    return [str(row["name"]) for row in cfg.get("selected_cases", [])]


def shell_centroids(nodes: dict[int, tuple[float, float, float]],
                    shells: dict[int, tuple[int, ...]]) -> dict[int, tuple[float, float]]:
    out = {}
    for eid, nids in shells.items():
        pts = [nodes[n] for n in nids if n in nodes]
        if not pts:
            continue
        x = sum(p[0] for p in pts) / len(pts)
        y = sum(p[1] for p in pts) / len(pts)
        out[eid] = (x, y)
    return out


def save_scatter(path: Path, xs: list[float], ys: list[float], values: list[float],
                 title: str, label: str) -> None:
    if not values:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 3.2), dpi=180)
    sc = ax.scatter(xs, ys, c=values, s=7, cmap="turbo")
    ax.set_title(title)
    ax.set_xlabel("x / mm")
    ax.set_ylabel("y / mm")
    ax.set_aspect("equal", adjustable="box")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(label)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def summarize_case(root: Path, output_root: Path, cfg: dict, case: str) -> list[dict[str, str]]:
    src = cfg["source"]
    out = cfg["output"]
    model = root / src["model_root"] / f"v42_{case}.k"
    nodes, shells = parse_keyword_mesh(model)
    centroids = shell_centroids(nodes, shells)
    field_dir = output_root / out["field_root"] / case
    figure_dir = output_root / out["figure_root"] / case
    rows: list[dict[str, str]] = []

    for state_dir in sorted(field_dir.glob("state_*")):
        if not state_dir.is_dir():
            continue
        state = int(state_dir.name.split("_")[-1])
        disp = parse_nodal_displacement(state_dir / "displacement.k")
        if disp:
            ids = list(disp)
            magnitudes = [math.sqrt(dx * dx + dy * dy + dz * dz) for dx, dy, dz in disp.values()]
            max_i = max(range(len(ids)), key=lambda i: magnitudes[i])
            xs = [nodes[nid][0] for nid in ids if nid in nodes]
            ys = [nodes[nid][1] for nid in ids if nid in nodes]
            vals = [magnitudes[i] * 1.0e6 for i, nid in enumerate(ids) if nid in nodes]
            save_scatter(
                figure_dir / f"state_{state:04d}_result_displacement_nm.png",
                xs, ys, vals,
                f"{case} state {state}: real resultant displacement",
                "resultant displacement / nm",
            )
            rows.append({
                "case": case,
                "state": str(state),
                "field": "resultant_displacement",
                "result_label": "NODAL_DISPLACEMENT",
                "min_value": f"{min(magnitudes):.8e}",
                "max_value": f"{max(magnitudes):.8e}",
                "unit": "mm",
                "max_id": str(ids[max_i]),
                "figure": str(figure_dir / f"state_{state:04d}_result_displacement_nm.png"),
            })

        for scalar_file in sorted(state_dir.glob("*.k")):
            if scalar_file.name == "displacement.k":
                continue
            label, values = parse_scalar_results(scalar_file)
            if not values:
                continue
            ids = list(values)
            vals = [values[i] for i in ids]
            max_i = max(range(len(ids)), key=lambda i: vals[i])
            xs = [centroids[eid][0] for eid in ids if eid in centroids]
            ys = [centroids[eid][1] for eid in ids if eid in centroids]
            plot_vals = [values[eid] for eid in ids if eid in centroids]
            save_scatter(
                figure_dir / f"state_{state:04d}_{scalar_file.stem}.png",
                xs, ys, plot_vals,
                f"{case} state {state}: real {label}",
                label,
            )
            rows.append({
                "case": case,
                "state": str(state),
                "field": scalar_file.stem,
                "result_label": label,
                "min_value": f"{min(vals):.8e}",
                "max_value": f"{max(vals):.8e}",
                "unit": "as exported by LS-PrePost",
                "max_id": str(ids[max_i]),
                "figure": str(figure_dir / f"state_{state:04d}_{scalar_file.stem}.png"),
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
    requested = case_names(cfg) if args.case == "all" else [args.case]

    rows: list[dict[str, str]] = []
    for case in requested:
        rows.extend(summarize_case(root, output_root, cfg, case))

    out_dir = output_root / cfg["output"]["lightweight_root"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "v44_real_field_summary.csv"
    readme = out_dir / "README.md"

    fieldnames = ["case", "state", "field", "result_label", "min_value", "max_value", "unit", "max_id", "figure"]
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        wr.writerows(rows)

    readme.write_text(
        "# V44 real field extraction summary\n\n"
        "This folder summarizes text field files exported from real d3plot data by LS-PrePost.\n\n"
        f"- exported rows: {len(rows)}\n"
        "- summary CSV: `v44_real_field_summary.csv`\n"
        "- full field files stay in `field_data/v44_lsprepost_real_field_extract/`\n"
        "- generated field figures stay in `figures/v44_lsprepost_real_field_extract/`\n",
        encoding="utf-8",
    )
    print(f"[OK] {out_csv}")
    print(f"[OK] {readme}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

