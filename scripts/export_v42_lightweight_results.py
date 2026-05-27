"""
export_v42_lightweight_results.py
=================================

Export small V4C recoil-sweep summaries from ignored LS-DYNA results.

Run from project root:
    python scripts\\export_v42_lightweight_results.py
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def parse_messag(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {
            "status": "missing messag",
            "normal_termination": "no",
            "problem_time_ms": "",
            "problem_cycle": "",
        }
    raw = path.read_text(encoding="utf-8", errors="ignore")
    collapsed = re.sub(r"\s+", "", raw).lower()
    normal = "normaltermination" in collapsed
    error = "errortermination" in collapsed
    times = re.findall(r"Problem time\s*=\s*([0-9.Ee+-]+)", raw)
    cycles = re.findall(r"Problem cycle\s*=\s*([0-9]+)", raw)
    return {
        "status": "OK" if normal else ("error termination" if error else "not completed"),
        "normal_termination": "yes" if normal else "no",
        "problem_time_ms": times[-1] if times else "",
        "problem_cycle": cycles[-1] if cycles else "",
    }


def d3plot_stats(result_dir: Path) -> dict[str, str]:
    files = sorted(result_dir.glob("d3plot*"))
    total_bytes = sum(p.stat().st_size for p in files if p.is_file())
    return {
        "d3plot_files": str(len(files)),
        "d3plot_total_MB": f"{total_bytes / 1_000_000:.3f}",
    }


def main() -> int:
    default_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", type=Path, default=default_root)
    ap.add_argument("--output-root", type=Path, default=None)
    args = ap.parse_args()
    root = args.project_root.resolve()
    output_root = args.output_root.resolve() if args.output_root else root

    registry = root / "models" / "v42_dynamic_recoil_sweep" / "v42_case_registry.csv"
    if not registry.is_file():
        raise SystemExit(f"registry not found: {registry}")

    out_dir = output_root / "lightweight_results" / "v42_dynamic_recoil_sweep"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "v42_case_summary.csv"
    out_md = out_dir / "README.md"

    rows = []
    for reg in read_csv_rows(registry):
        case = reg["name"]
        result_dir = root / "results" / "v42_dynamic_recoil_sweep" / case
        msg = parse_messag(result_dir / "messag")
        plots = d3plot_stats(result_dir)
        rows.append({
            "case": case,
            "Ep_uJ": reg.get("Ep_uJ", ""),
            "recoil_velocity_m_s": reg.get("recoil_velocity_m_s", ""),
            "vapor_radius_um": reg.get("vapor_radius_um", ""),
            "removed_elements": reg.get("removed_elements", ""),
            "recoil_nodes": reg.get("recoil_nodes", ""),
            **msg,
            **plots,
            "notes": reg.get("notes", ""),
        })

    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    ok = sum(1 for r in rows if r["normal_termination"] == "yes")
    out_md.write_text(
        "# V4C recoil-sweep lightweight results\n\n"
        "Small summary exported from ignored `results/v42_dynamic_recoil_sweep/` outputs.\n\n"
        f"- cases in registry: {len(rows)}\n"
        f"- normal termination: {ok}\n"
        f"- summary CSV: `v42_case_summary.csv`\n\n"
        "Large LS-DYNA outputs remain in `results/` and should not be committed.\n",
        encoding="utf-8",
    )
    print(f"[OK] {out_csv}")
    print(f"[OK] {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
