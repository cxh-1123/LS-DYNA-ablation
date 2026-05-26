# Lightweight Results

This directory stores small, reviewable summaries exported from ignored
`results/` folders. It should contain CSV/JSON/Markdown only, not large
LS-DYNA solver outputs.

Regenerate after running V17/V26 post-processing:

```powershell
python scripts\export_lightweight_results.py
```

Expected exported files after real V17/V26 runs:

- `lightweight_results/v17_30ps_local/v17_case_summary.csv`
- `lightweight_results/v26_30ps_threshold_ablation/v26_case_summary.csv`
- `lightweight_results/v26_30ps_threshold_ablation/<case>/v26_threshold_metrics.csv`
- `lightweight_results/manifest.json`

Do not commit raw `d3plot`, `d3hsp`, `messag`, or `tprint` files.
