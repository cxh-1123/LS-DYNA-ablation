# V5A explicit ejecta marker plume pilot

V5A is the first step from "solid-only recoil response" toward a visible
ablation plume.

It keeps the V4C idea:

- pre-deleted vapor-zone crater
- recoil velocity on the exposed crater boundary
- real d3plot pressure/displacement output

It adds:

- a separate ejecta marker part above the crater
- small explicit Lagrangian shell patches as visible plume markers
- upward and radial initial velocity profiles

This is intentionally a robust pilot, not a final SPH/ALE gas model.  The
marker patches are there to make the ejection path visible in d3plot and to
build the post-processing/figure workflow.  After this runs cleanly, the marker
patches can be replaced with true SPH particles or an ALE gas/plume region.

## Build

```powershell
python scripts\build_v50_dynamic_ejecta_particle_pilot.py
```

## Dry run

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v50_selected.ps1 -DryRun
```

## Run

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v50_selected.ps1 -OnlyCase ep_5p0uj_ejecta_8500ms
```

## View

Open:

```text
results\v50_dynamic_ejecta_particle_pilot\ep_5p0uj_ejecta_8500ms\d3plot
```

In LS-PrePost, turn on all parts and animate.  The ejecta marker part should
rise above the crater while the substrate still shows pressure and displacement
response.

## Paper-style reference figure (V5B unified sequence)

Unified early ablation + plume diffusion on **one fixed coordinate frame**
(r=[-60,60] um, z=[-8,80] um for main sequence).  Not stitched plume/substrate
panels.

```powershell
python scripts\build_v3b_plume_model.py
python scripts\plot_v5b_unified_ablation_sequence.py
```

Outputs:

```text
figures\v5b_reference_style\ep_5p0uj_ejecta_8500ms\
  v5b_early_ablation_plume_sequence.png   # 0,10,30,50,100,200,500 ps (2x4)
  v5b_late_ablation_plume_sequence.png    # 1000,2000,5000 ps supplement
  v5b_crater_evolution_summary.png        # crater vs time (separate)
```

| Layer | Role |
| --- | --- |
| V2.6 | Crater radius/depth vs time (surface bowl) |
| V3B.1 | Coupled plume density from crater mouth |
| V5A | Real d3plot markers (LS-PrePost only, not in main sequence) |

## Roadmap toward full reference physics

| Stage | Content |
| --- | --- |
| **V5A** (current) | Explicit ejecta marker shells + recoil + pre-deleted crater |
| **V5B** (this) | Reference-style composite visualization (V3B + V5A) |
| **V5C** | Denser markers / massless tracer particles in LS-DYNA |
| **V6** | SPH or ALE gas plume (true multiphase, TTM-coupled if needed) |

## Rebuild after f10 fix

If LS-DYNA reported `*INITIAL_VELOCITY_NODE` format errors (`0.0-5.000e+03`),
rebuild and re-run:

```powershell
python scripts\build_v50_dynamic_ejecta_particle_pilot.py
powershell -ExecutionPolicy Bypass -File scripts\run_v50_selected.ps1 -OnlyCase ep_5p0uj_ejecta_8500ms
```

