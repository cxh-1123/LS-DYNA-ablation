# V6A LAMMPS TTM-MD Pilot

This stage starts the move from LS-DYNA visual/proxy ejecta models to atomistic
LAMMPS TTM-MD.

## What V6A Is

V6A is a small silicon slab test.  It is designed to answer four basic
questions before running expensive ablation cases:

- Can the local LAMMPS installation run a silicon MD slab?
- Can it find the Stillinger-Weber Si potential file?
- Was LAMMPS built with the `EXTRA-FIX` package needed for `fix ttm`?
- Can the workflow write atom snapshots and electron-temperature grid files?

## Why This Is Not Yet the Final Paper Figure

The current V5A LS-DYNA ejecta-marker result is a visualization proxy.  V6A is
the first physics-based atomistic replacement.  It does not yet try to match the
full paper image directly.  It first verifies the smallest reliable TTM-MD unit:
atoms plus an electronic heat bath.

## Physical Meaning

The model uses:

- MD atoms for the silicon lattice;
- a TTM electron grid for electronic temperature;
- electron-lattice coupling through LAMMPS `fix ttm`;
- a hot near-surface electronic initial condition for the first pilot laser
  heating check.

LAMMPS `fix ttm` requires a 3D periodic orthogonal simulation box.  Therefore
V6A uses a slab plus a vacuum gap inside a periodic box.  For production laser
ablation, V6B should move to `fix ttm/mod` or `fix ttm/thermal` after confirming
the local LAMMPS version supports it.

## Files

- `config/v60_lammps_ttm_md_pilot.toml`: model and case settings.
- `scripts/build_v60_lammps_ttm_md_pilot.py`: writes LAMMPS input files.
- `scripts/run_v60_selected.ps1`: launches selected LAMMPS cases.
- `models/v60_lammps_ttm_md_pilot/v60_case_registry.csv`: run registry.
- `models/v60_lammps_ttm_md_pilot/v60_ttm_mod_parameters.txt`: V6B starting
  parameter file for `fix ttm/mod` calibration.

## First Run Order

Build:

```powershell
python scripts\build_v60_lammps_ttm_md_pilot.py
```

Dry check:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v60_selected.ps1 -DryRun
```

Run only the non-laser MD check first:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v60_selected.ps1 -OnlyCase equil_300K
```

If that finishes, run the TTM pilot:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v60_selected.ps1 -OnlyCase ttm_initial_pulse_5uj_equiv
```

If LAMMPS is not on PATH, set:

```powershell
$env:LAMMPS_EXE="C:\path\to\lmp.exe"
```

## References

- LAMMPS `fix ttm` / `fix ttm/mod` documentation:
  https://docs.lammps.org/fix_ttm.html
- LAMMPS Stillinger-Weber potential documentation:
  https://docs.lammps.org/pair_sw.html
- LAMMPS dump documentation:
  https://docs.lammps.org/dump.html
- D. M. Duffy and A. M. Rutherford, J. Phys.: Condens. Matter 19, 016207
  (2007): https://doi.org/10.1088/0953-8984/19/1/016207

