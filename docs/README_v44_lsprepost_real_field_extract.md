# V44 real d3plot field extraction with LS-PrePost

V44 connects the current V4C d3plot results to repeatable post-processing.

The goal is to stop relying on screenshots or proxy plots.  The preferred path
uses `lasso-python` to read real `d3plot` arrays directly.  A secondary
LS-PrePost command-file path is kept for display automation, but LS-PrePost
4.12 may not write numeric field values from `output ...` in `-nographics`
mode.

## What it extracts

- `displacement`: real `node_displacement` vectors from d3plot.
- `velocity`: real `node_velocity` vectors from d3plot, when present.
- `pressure`: computed from real shell stress as `-(sxx + syy + szz) / 3`.

The V4C dynamic decks do not contain a coupled temperature history or real
plume particles.  Temperature must be taken from the thermal V1.7/V2.6 outputs,
or added later through a coupled thermal-dynamic model.  Plume needs SPH/ALE or
a particle proxy model.

## Run

Install the parser once:

```powershell
python -m pip install lasso-python
```

Then run from the repository root:

```powershell
python scripts\extract_v44_lasso_d3plot.py --case ep_5p0uj_recoil_5000ms
```

For all configured V4C cases:

```powershell
python scripts\extract_v44_lasso_d3plot.py
```

## LS-PrePost display batch

This path can still generate LS-PrePost command files for repeatable viewing:

```powershell
$env:LSPREPOST_EXE="C:\path\to\lsprepost.exe"
powershell -ExecutionPolicy Bypass -File scripts\run_v44_lsprepost_export.ps1 -Case ep_5p0uj_recoil_5000ms
```

## Outputs

- command files: `post/v44_lsprepost_real_field_extract/`
- raw exported fields: `field_data/v44_lsprepost_real_field_extract/`
- CSV summary: `lightweight_results/v44_lsprepost_real_field_extract/v44_real_field_summary.csv`
- field figures: `figures/v44_lsprepost_real_field_extract/`

## If pressure exports the wrong field

Different LS-PrePost builds can use different numeric fringe ids.  If the
exported pressure file is not labeled as pressure, do this once:

1. Open one V4C `d3plot` manually.
2. Select `Fcomp -> Misc -> pressure -> Apply`.
3. Exit LS-PrePost and find the generated `lspost.cfile`.
4. Copy the recorded `fringe ...` line into
   `config/v44_lsprepost_real_field_extract.toml` under `pressure.fringe_command`.
5. Re-run V44.

The nodal displacement export does not need this fringe calibration.
