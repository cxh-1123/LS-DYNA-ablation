# V4.0A dynamic ablation pilot

V4.0A is the first step from threshold post-processing toward dynamic ablation.
It is not the final calibrated ablation model yet.

## What changes from V1.7/V2.6

Earlier stages:

- V1.7 solves temperature.
- V2.6 marks regions where `T >= T_vap`.
- V3B shows a proxy plume/shock visualization.

V4.0A starts a dynamic LS-DYNA workflow:

- reuse the V1.7 local refined mesh;
- switch to structural dynamics material cards;
- estimate the vapor-zone radius from the V1.7 vapor threshold;
- assign a first ejecta impulse to the near-surface vapor zone;
- export candidate elements for active material deletion.

Plain-language meaning:

> V2.6 says "this area is hot enough to disappear."  V4.0A begins turning that
> into a dynamic model where the hot surface zone can move and later be deleted.

## Current pilot assumption

The V1.7 vapor threshold estimate is:

```text
Ep_vap ~= 3.052 uJ
```

For a Gaussian spot, V4.0A estimates the vapor-zone radius by:

```text
r_vap = w0 / sqrt(2) * sqrt(log(Ep / Ep_vap))
```

This gives no vapor zone below threshold, a very small zone near `3.1 uJ`, and
a larger zone at `4-5 uJ`.

## Generated files

Run:

```powershell
python scripts\build_v40_dynamic_ablation_pilot.py
```

This writes:

```text
models/v40_dynamic_ablation_pilot/v40_<case>.k
models/v40_dynamic_ablation_pilot/v40_case_registry.csv
models/v40_dynamic_ablation_pilot/v40_erosion_candidates.csv
```

The `.k` files are structural dynamics pilot decks.  The CSV files record the
estimated vapor radius and the elements in the first near-surface layer that
would be candidates for active deletion.

## Important limitation

The generated V4.0A deck does not yet enable final active erosion.  That is
intentional.  In LS-DYNA, temperature-driven erosion is coupled to material
model choice and keyword-version details.  The safe workflow is:

1. generate this pilot deck;
2. run one case to check dynamic stability;
3. inspect motion/stress/velocity output;
4. enable active erosion with the exact keyword supported by the local
   LS-DYNA version;
5. rerun and compare crater growth with the V2.6 threshold estimate.

## Recommended first cases

Start with:

```powershell
.\scripts\run_v40_selected.ps1 -OnlyCase ep_3p1uj -Ncpu 4
```

Then run:

```powershell
.\scripts\run_v40_selected.ps1 -OnlyCase ep_4p0uj -Ncpu 4
```

`ep_3p1uj` checks the just-above-threshold behavior.  `ep_4p0uj` gives a more
visible dynamic response.

## What success looks like

For V4.0A, success means:

- LS-DYNA reaches normal termination;
- the top-center region moves upward rather than staying static;
- stress/velocity waves appear in the substrate;
- the candidate deletion region is small near `3.1 uJ` and larger at `4-5 uJ`;
- no large-file solver outputs are committed to GitHub.

Once this pilot is stable, the next step is V4.0B: active material deletion.
