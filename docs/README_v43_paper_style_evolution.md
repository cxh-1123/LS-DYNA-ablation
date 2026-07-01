# V4D paper-style evolution visualization

V4D generates paper-style figures for the current V4B/V4C ablation workflow.

## What it shows

The figures are designed to look closer to ultrafast ablation papers:

- time snapshots at `5 / 15 / 25 / 40 / 55 / 70 / 85 ps`;
- a pre-deleted crater at the silicon surface;
- proxy ejecta/plume particles colored by temperature;
- pressure/displacement trend panels from V4C recoil sweeps.

## Important limitation

V4D is a proxy visualization.  It does not yet read full d3plot binary fields.
It combines:

- V4C case metadata;
- V4C normal-termination summaries;
- LS-PrePost peak readings recorded in
  `config/v43_paper_style_evolution.toml`;
- a simple ballistic plume/ejecta visualization model.

Use it for communication, debugging, and planning the next physics step.  Do not
call it calibrated TTM-MD or fully coupled plume dynamics.

## Generate

```powershell
python scripts\plot_v43_paper_style_evolution.py
```

Outputs:

```text
figures/v43_paper_style_evolution/v43_snapshot_sequence.png
figures/v43_paper_style_evolution/v43_recoil_scaling_summary.png
figures/v43_paper_style_evolution/v43_proxy_frame_data.csv
```

## Next physics step

To make these figures more physical, the next stage should replace proxy fields
with one of:

- LS-PrePost batch exports from d3plot;
- a Python d3plot reader;
- LS-DYNA SPH/ALE ejecta particles;
- a coupled TTM-MD or TTM-continuum material-removal model.
