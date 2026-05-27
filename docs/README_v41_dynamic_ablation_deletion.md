# V4B dynamic ablation with material removal

V4B is the first material-removal model in this project.

## What V4B does

V4A kept the whole mesh and added an initial ejecta/recoil motion.  V4B goes one
step further: it removes the vapor-zone elements from the LS-DYNA structural
mesh before the calculation starts.

Plain-language meaning:

> V4A made the hot surface region move.  V4B actually takes the vaporized
> material out of the model, so the calculation starts with a small crater.

## Why pre-deletion first

Direct active erosion in LS-DYNA depends on the exact solver version, material
model, and keyword option.  A pre-deleted deck is more robust and easier to
debug:

- no version-specific erosion card is needed;
- the crater geometry is explicit in the input file;
- the remaining target can still respond dynamically;
- the result can be compared against V2.6 equivalent-crater dimensions.

Once V4B runs cleanly, V4C can replace pre-deletion with active erosion.

## Generated files

Run:

```powershell
python scripts\build_v41_dynamic_ablation_deletion.py
```

This writes:

```text
models/v41_dynamic_ablation_deletion/v41_<case>.k
models/v41_dynamic_ablation_deletion/v41_case_registry.csv
models/v41_dynamic_ablation_deletion/v41_removed_elements.csv
```

## First run

Start with the smallest deletion case:

```powershell
.\scripts\run_v41_selected.ps1 -OnlyCase ep_3p1uj -Ncpu 4
```

If it reaches normal termination, run:

```powershell
.\scripts\run_v41_selected.ps1 -OnlyCase ep_4p0uj -Ncpu 4
```

## What success means

For V4B, success means:

- LS-DYNA reads the deck without keyword errors;
- the model reaches normal termination;
- the top center elements listed in `v41_removed_elements.csv` are absent from
  the structural mesh;
- the remaining substrate responds dynamically around the crater.

Do not commit `results/` solver outputs to GitHub.  Commit only lightweight
CSV/Markdown summaries.
