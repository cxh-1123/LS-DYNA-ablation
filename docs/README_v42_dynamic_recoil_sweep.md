# V4C recoil-drive sweep

V4B confirmed that pre-deleted crater decks run normally, but the displacement
was very small.  For the `ep_5p0uj` case, the peak resultant displacement was
only about `5.3e-08 mm`, or roughly `0.053 nm`.

V4C keeps the same pre-deleted crater idea and scans stronger recoil velocities
on the newly exposed crater boundary.

## Why this step exists

The goal is to answer a practical question:

> Is the weak visible motion caused by the model form, or simply by too small a
> recoil/ejecta drive?

V4C does not yet claim calibrated vapor pressure.  It is a sensitivity sweep.

## Cases

The first sweep uses `ep_5p0uj` and the same vapor-threshold crater as V4B:

| Case | Recoil velocity |
| --- | ---: |
| `ep_5p0uj_recoil_100ms` | `-100 m/s` |
| `ep_5p0uj_recoil_1000ms` | `-1000 m/s` |
| `ep_5p0uj_recoil_3000ms` | `-3000 m/s` |
| `ep_5p0uj_recoil_5000ms` | `-5000 m/s` |

Negative y velocity points downward into the target in the current 2D r-z
section.  This represents a recoil-like impulse on the crater boundary.

## Generate

```powershell
python scripts\build_v42_dynamic_recoil_sweep.py
```

This writes:

```text
models/v42_dynamic_recoil_sweep/v42_<case>.k
models/v42_dynamic_recoil_sweep/v42_case_registry.csv
models/v42_dynamic_recoil_sweep/v42_removed_elements.csv
```

## Run

Start from the middle case:

```powershell
.\scripts\run_v42_selected.ps1 -OnlyCase ep_5p0uj_recoil_1000ms -Ncpu 4
```

If stable, run:

```powershell
.\scripts\run_v42_selected.ps1 -OnlyCase ep_5p0uj_recoil_3000ms -Ncpu 4
.\scripts\run_v42_selected.ps1 -OnlyCase ep_5p0uj_recoil_5000ms -Ncpu 4
```

## Interpret

Use `Ndv -> result displacement` and `Misc -> pressure` in LS-PrePost.

If displacement and pressure scale roughly with recoil velocity, then the weak
V4B motion was mostly a drive-strength issue.  If the response remains tiny,
the next step should be changing the loading form, for example adding a
distributed recoil pressure instead of only nodal initial velocity.
