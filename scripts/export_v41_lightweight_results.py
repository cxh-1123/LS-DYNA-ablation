"""
export_v41_lightweight_results.py
=================================

Export small V4B run summaries from ignored LS-DYNA results.

Run from project root:
    python scripts\\export_v41_lightweight_results.py

Or summarize another checkout:
    python scripts\\export_v41_lightweight_results.py --project-root D:\\cxh-daima\\LS-DYNA-ablation

Or read results from another checkout but write summaries here:
    python scripts\\export_v41_lightweight_results.py --project-root D:\\cxh-daima\\LS-DYNA-ablation --output-root .
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
    if normal:
        status = "OK"
    elif error:
        status = "error termination"
    else:
        status = "not completed"
    return {
        "status": status,
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

    registry = root / "models" / "v41_dynamic_ablation_deletion" / "v41_case_registry.csv"
    if not registry.is_file():
        raise SystemExit(f"registry not found: {registry}")

    out_dir = output_root / "lightweight_results" / "v41_dynamic_ablation_deletion"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "v41_case_summary.csv"
    out_md = out_dir / "README.md"

    rows = []
    for reg in read_csv_rows(registry):
        case = reg["name"]
        result_dir = root / "results" / "v41_dynamic_ablation_deletion" / case
        msg = parse_messag(result_dir / "messag")
        plots = d3plot_stats(result_dir)
        rows.append({
            "case": case,
            "Ep_uJ": reg.get("Ep_uJ", ""),
            "vapor_radius_um": reg.get("vapor_radius_um", ""),
            "removed_elements": reg.get("removed_elements", ""),
            "kept_elements": reg.get("kept_elements", ""),
            "recoil_nodes": reg.get("recoil_nodes", ""),
            **msg,
            **plots,
            "notes": reg.get("notes", ""),
        })

    fieldnames = list(rows[0].keys())
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        wr.writerows(rows)

    ok = sum(1 for r in rows if r["normal_termination"] == "yes")
    out_md.write_text(
        "# V4B lightweight results\n\n"
        "Small summary exported from ignored `results/v41_dynamic_ablation_deletion/` outputs.\n\n"
        f"- cases in registry: {len(rows)}\n"
        f"- normal termination: {ok}\n"
        f"- summary CSV: `v41_case_summary.csv`\n\n"
        "Large LS-DYNA outputs such as `d3plot*`, `messag`, `glstat`, and `matsum` "
        "remain in `results/` and should not be committed.\n",
        encoding="utf-8",
    )

    print(f"[OK] {out_csv}")
    print(f"[OK] {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
