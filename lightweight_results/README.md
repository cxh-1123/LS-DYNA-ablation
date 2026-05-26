# Lightweight Results

This directory stores small, reviewable summaries exported from ignored
`results/` folders. It should contain CSV/JSON/Markdown only, not large
LS-DYNA solver outputs.

Regenerate after running V17/V26 post-processing:

```powershell
python scripts\export_lightweight_results.py
```

Exported files:
- `lightweight_results/v17_30ps_local/v17_case_summary.csv` (1488 bytes, 11 rows)
- `lightweight_results/v26_30ps_threshold_ablation/v26_case_summary.csv` (900 bytes, 7 rows)
- `lightweight_results/v26_30ps_threshold_ablation/ep_1p0uj/v26_threshold_metrics.csv` (12520 bytes, 177 rows)
- `lightweight_results/v26_30ps_threshold_ablation/ep_1p5uj/v26_threshold_metrics.csv` (12857 bytes, 178 rows)
- `lightweight_results/v26_30ps_threshold_ablation/ep_2p0uj/v26_threshold_metrics.csv` (13006 bytes, 180 rows)
- `lightweight_results/v26_30ps_threshold_ablation/ep_2p5uj/v26_threshold_metrics.csv` (13051 bytes, 181 rows)
- `lightweight_results/v26_30ps_threshold_ablation/ep_3p0uj/v26_threshold_metrics.csv` (13340 bytes, 183 rows)
- `lightweight_results/v26_30ps_threshold_ablation/ep_4p0uj/v26_threshold_metrics.csv` (14157 bytes, 186 rows)
- `lightweight_results/v26_30ps_threshold_ablation/ep_5p0uj/v26_threshold_metrics.csv` (14891 bytes, 188 rows)
