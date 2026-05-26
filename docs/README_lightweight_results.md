# Lightweight result package

Large LS-DYNA outputs stay under `results/` and are ignored by Git.  Small
review files are exported to `lightweight_results/`:

- `v17_30ps_local/v17_case_summary.csv`
- `v26_30ps_threshold_ablation/v26_case_summary.csv`
- `v26_30ps_threshold_ablation/<case>/v26_threshold_metrics.csv`
- `manifest.json`

Regenerate the package after V17/V26 checks:

```powershell
python scripts\check_v17_outputs.py
python scripts\extract_v26_30ps_threshold_ablation.py
python scripts\export_lightweight_results.py
```

Commit only `lightweight_results/` summaries, not `results/` solver files.
