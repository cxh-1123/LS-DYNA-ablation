# Lightweight Results

This directory stores small, reviewable summaries exported from ignored
`results/` folders. It should contain CSV/JSON/Markdown only, not large
LS-DYNA solver outputs.

Regenerate after running V17/V26 post-processing:

```powershell
python scripts\export_lightweight_results.py
```

Exported files:
- `lightweight_results/v17_30ps_local/v17_case_summary.csv` (1498 bytes, 11 rows)
- `lightweight_results/v26_30ps_threshold_ablation/v26_case_summary.csv` (541 bytes, 4 rows)

Missing at export time:
- `results/v26_30ps_threshold_ablation/*/v26_threshold_metrics.csv`
