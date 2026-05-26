"""
Export small, reviewable result summaries from ignored LS-DYNA result folders.

The project intentionally ignores ``results/`` because it contains large
regenerable solver outputs.  This script copies only lightweight CSV / JSON /
Markdown summaries into ``lightweight_results/`` so GitHub can keep the key
numbers needed for review.

Run from project root:
    python scripts/export_lightweight_results.py
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_ROOT = PROJECT_ROOT / "results"
OUT_ROOT = PROJECT_ROOT / "lightweight_results"


SUMMARY_TARGETS = [
    RESULTS_ROOT / "v17_30ps_local" / "v17_case_summary.csv",
    RESULTS_ROOT / "v26_30ps_threshold_ablation" / "v26_case_summary.csv",
]

METRICS_GLOBS = [
    "v26_30ps_threshold_ablation/*/v26_threshold_metrics.csv",
]


@dataclass
class ExportedFile:
    source: str
    destination: str
    bytes: int
    rows: int | None


def count_csv_rows(path: Path) -> int | None:
    try:
        with path.open("r", encoding="utf-8", newline="") as fh:
            return max(sum(1 for _ in csv.reader(fh)) - 1, 0)
    except Exception:
        return None


def relative(path: Path) -> str:
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def copy_file(src: Path, dst: Path) -> ExportedFile:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return ExportedFile(
        source=relative(src),
        destination=relative(dst),
        bytes=dst.stat().st_size,
        rows=count_csv_rows(dst) if dst.suffix.lower() == ".csv" else None,
    )


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    exported: list[ExportedFile] = []
    missing: list[str] = []

    for src in SUMMARY_TARGETS:
        if src.is_file():
            exported.append(copy_file(src, OUT_ROOT / relative(src).removeprefix("results/")))
        else:
            missing.append(relative(src))

    for pattern in METRICS_GLOBS:
        matches = sorted(RESULTS_ROOT.glob(pattern))
        if not matches:
            missing.append(f"results/{pattern}")
        for src in matches:
            exported.append(copy_file(src, OUT_ROOT / relative(src).removeprefix("results/")))

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Small reviewable summaries copied from ignored results/ outputs.",
        "exported_files": [asdict(item) for item in exported],
        "missing_expected_files": missing,
        "note": (
            "Large LS-DYNA outputs such as d3plot, d3hsp, messag and tprint "
            "remain ignored and are not part of this package."
        ),
    }
    (OUT_ROOT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    readme_lines = [
        "# Lightweight Results",
        "",
        "This directory stores small, reviewable summaries exported from ignored",
        "`results/` folders. It should contain CSV/JSON/Markdown only, not large",
        "LS-DYNA solver outputs.",
        "",
        "Regenerate after running V17/V26 post-processing:",
        "",
        "```powershell",
        "python scripts\\export_lightweight_results.py",
        "```",
        "",
        "Exported files:",
    ]
    if exported:
        readme_lines.extend(
            f"- `{item.destination}` ({item.bytes} bytes"
            + (f", {item.rows} rows" if item.rows is not None else "")
            + ")"
            for item in exported
        )
    else:
        readme_lines.append("- None yet. Run the V17/V26 checks first.")

    if missing:
        readme_lines.extend(["", "Missing at export time:"])
        readme_lines.extend(f"- `{item}`" for item in missing)

    (OUT_ROOT / "README.md").write_text("\n".join(readme_lines) + "\n", encoding="utf-8")

    print("=" * 72)
    print("Lightweight result export")
    print(f"  output root : {OUT_ROOT}")
    print(f"  exported    : {len(exported)} file(s)")
    print(f"  missing     : {len(missing)} expected path(s)")
    print("=" * 72)
    for item in exported:
        row_note = f", rows={item.rows}" if item.rows is not None else ""
        print(f"  [OK] {item.destination} ({item.bytes} bytes{row_note})")
    for item in missing:
        print(f"  [MISS] {item}")
    return 0 if exported else 1


if __name__ == "__main__":
    raise SystemExit(main())
