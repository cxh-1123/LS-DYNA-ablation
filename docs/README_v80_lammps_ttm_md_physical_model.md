# V8 Si TTM-MD Physical Model Update

V8 is the first step after the V7 continuous-field diagnostics.  The goal is
not to make prettier plots from a weak trajectory, but to correct the model and
diagnostic pipeline so the Si ultrafast-ablation results are physically
interpretable.

## V7 Findings

The completed `fluence_0p5_100ps` trajectory mainly showed rapid heating,
lattice disordering, and a few high-energy atoms above the surface.  It did not
yet show a stable plume.

Key warnings from V7:

- Raw maximum-z lift reached several nm, but robust p95 surface lift stayed near
  zero.  The raw lift was dominated by a few outlier atoms.
- Many atoms exceeded the melt/vapor temperature proxy, but detached and stable
  ejecta counts remained small.
- The largest connected component still contained almost all atoms, so the bulk
  remained connected.
- The old density field was based too heavily on neighbor counts and could show
  unphysical high-density regions during expansion.
- The old pressure field was a proxy and could not be trusted as a direct GPa
  result.
- The old atom temperature proxy mixed thermal motion with drift/expansion
  motion.
- The V6B debug cell was too narrow laterally and the 100 ps window was too
  short for stable ejecta/plume formation.

## V8 Corrections

### Density

The main continuous-field density is now `bin_density`:

`rho_bin/rho0 = current atom count in an x-z bin / initial atom count in the same bin`

Neighbor density is kept only for connectivity and ejecta decisions:

`neighbor_density = neighbor_count / initial_neighbor_count`

This avoids using a neighbor-count field as the main density map.  Initial empty
or low-count bins are masked.  The diagnostics now report:

- `low_density_area_fraction`, where `rho_bin/rho0 < 0.7`
- `void_area_fraction`, where `rho_bin/rho0 < 0.3`
- initial density sanity check, expected near 1

### Atom Temperature

The V8 postprocessor computes both raw and thermal atom temperatures.  Thermal
temperature subtracts the local mean velocity before converting kinetic energy
to temperature:

`v_thermal = v_atom - mean(v_local_bin)`

The main physical plots should use `atom_T_thermal` or TTM lattice temperature
instead of the raw kinetic-temperature proxy.

### Pressure And Stress

LAMMPS `units metal` reports `compute stress/atom` as pressure times volume.
V8 divides binned virial stress by the x-z control volume and converts bar to
GPa.  Low-count bins are masked, and p1/p99 clipping is used for visualization.

This is more physical than V7, but it should still be treated as a
coarse-grained stress estimate.  If the bin volume or stress output is changed,
re-check the conversion before using GPa in final figures.

### Surface And Crater Metrics

V8 keeps the robust surface definition from V7:

- raw maximum-z lift is retained as `raw_lift_outlier_sensitive`
- the main surface metric is `robust_p95_lift`
- crater depth uses the robust p95 surface after excluding ejecta candidates
- `surface_profile_evolution.png` is generated for direct inspection

If raw lift is large while robust lift remains near zero, that means a few atoms
escaped but the surface as a whole did not lift or crater significantly.

### Ejecta Classification

Stable ejecta is no longer temperature-only.  It must satisfy height, upward
velocity, low density or disconnection from the bulk, and persistence for at
least three frames.

Reported classes:

- upward atoms
- vapor candidates
- detached candidates
- stable ejecta
- cluster ejecta
- largest ejecta cluster size
- ejected mass fraction

## Laser Source

V8 distinguishes incident fluence and absorbed fluence:

`absorbed_fluence = incident_fluence * (1 - R)`

Current baseline:

- incident fluence: configured per case
- reflectivity `R = 0.35`
- absorbed fraction: `0.65`
- optical penetration depth: `10 nm`
- temporal pulse: Gaussian, 30 ps FWHM
- pulse center: 15 ps

The 100 ps or 300 ps values are observation windows, not laser pulse durations.

Generated laser checks:

- `figures/v80_physical_model_update/laser_source/laser_temporal_profile.png`
- `figures/v80_physical_model_update/laser_source/laser_source_z_profile.png`
- `figures/v80_physical_model_update/laser_source/absorbed_energy_check.csv`

For the free-surface ablation case, `fix ttm/mod` is not used.  Instead, V8
uses an equivalent Beer-Lambert lattice heat source applied to near-surface
layers.  This is a controlled approximation for morphology/ejecta testing, not
a full two-temperature electron transport solution.  The deposited lattice
fraction is written explicitly as
`free_surface_lattice_deposition_fraction`.

The free-surface debug cases do not use inline clipping or percentile variables
inside the LAMMPS input.  This keeps the input compatible with the currently
installed LAMMPS build, which rejected an inline clipping expression in an
equal-style variable at step 0.  If the smooth source still produces max-KE or
max-vz spikes, reduce
`free_surface_lattice_deposition_fraction` or split the heat deposition into a
segmented run rather than adding unsupported variable formulas.

## TTM And Potential Warnings

The current TTM parameter block follows the LAMMPS `Si.ttm_mod` baseline.  It is
a transparent starting point, not a final literature-validated parameter set.
Before final paper use, replace or validate:

- electron heat capacity model
- electron diffusivity / conductivity
- electron-phonon coupling
- optical penetration depth for the actual laser wavelength
- reflectivity for the actual sample and wavelength

The current Si potential is Stillinger-Weber.  It is useful for crystalline and
liquid Si baseline behavior, but high-temperature vapor, ionization, nonthermal
bond softening, and cluster emission are not guaranteed to be reliable.  V8
therefore keeps `enable_nonthermal_softening = false` by default and documents
it as a phenomenological option only.

## System Sizes

V8 adds three size presets:

- `debug`: fast checks only; not allowed as a paper main result
- `small_physical`: about 10 nm lateral size, 30 nm Si depth, and 24 nm vacuum
- `paper_style`: more expensive candidate for paper-style results

The first real run should use `small_physical`.  It is much larger than V6B and
is intended to allow more realistic pressure-wave propagation and early ejecta
formation.

## Boundary Condition And Case Separation

LAMMPS `fix ttm/mod` requires periodic simulation-box boundaries.  A direct
`p p f` free-surface box will stop with:

`Cannot use non-periodic boundaries with fix ttm/mod`

Therefore V8 has two explicitly separated cases:

1. `fluence_0p5_30ps_ttm_periodic_slab_diagnostic`
   - boundary: `p p p`
   - uses `fix ttm/mod`
   - keeps a top vacuum gap, but the box is periodic
   - runs 30 ps only
   - diagnostic only: verifies electron-lattice coupling and temperature
     evolution
   - not used for strict free-surface crater, ejecta, plume or morphology
     conclusions

2. `fluence_0p5_300ps_free_surface_ablation`
   - boundary: `p p fs`
   - does not use `fix ttm/mod`
   - uses equivalent Beer-Lambert near-surface lattice heating
   - runs 300 ps
   - main result for surface morphology, crater, ejecta marker, low-density
     plume region and ablation overview

An intermediate free-surface check is also generated:

3. `fluence_0p5_20ps_free_surface_ablation_debug`
   - boundary: `p p fs`
   - does not use `fix ttm/mod`
   - timestep: `0.0005 ps`
   - final time: `20 ps`
   - short stability test after the 15 ps lost-atoms failure
   - primary checks: no lost atoms, no sudden max-KE/max-vz spike, no pressure
     or temperature explosion

4. `fluence_0p5_50ps_free_surface_ablation_stable_debug`
   - boundary: `p p fs`
   - does not use `fix ttm/mod`
   - timestep: `0.0005 ps`
   - final time: `50 ps`
   - run only if the 20 ps case is stable

5. `fluence_0p5_100ps_free_surface_ablation_quickcheck`
   - boundary: `p p fs`
   - does not use `fix ttm/mod`
   - timestep: `0.0005 ps` in the current stability-first setup
   - uses the same equivalent Beer-Lambert near-surface lattice heating as the
     300 ps case
   - runs 100 ps
   - intended to check free-surface motion, early crater trend, ejecta
     classification and the V8 postprocessing chain before spending time on
     the 300 ps production run

For the free-surface cases, the lower z face is fixed and the upper z face is
shrink-wrapped.  This keeps the substrate boundary controlled while preventing
early ejecta atoms from being deleted as `Lost atoms` at the top of the box.
This boundary setting does not by itself prove the physics is stable: if max
atom KE, max upward velocity or pressure jumps sharply, the source/timestep
model is still unstable and must be corrected before longer runs.

Do not present the periodic-slab TTM diagnostic as a real free-surface ablation
result.

## Time Windows

V8 uses two time windows:

- TTM diagnostic: 30 ps
- free-surface ablation morphology case: 300 ps

- 0-30 ps: laser energy deposition into the electron system
- 30-100 ps: lattice heating, melting/disordering, pressure release
- 100-300 ps: possible ejecta growth, plume onset, and surface relaxation

The generated LAMMPS input resets time to zero before the production stage.
Atom dumps are split into:

- early dump: 0-100 ps, high-frequency output
- late dump: 100-300 ps, lower-frequency output

## First V8 Run Order

Run in this order:

1. `fluence_0p5_30ps_ttm_periodic_slab_diagnostic`
2. `fluence_0p5_20ps_free_surface_ablation_debug`
3. postprocess the 20 ps debug case and inspect `stability_diagnostics_v8.png`
   plus `heat_group_diagnostics_v8.csv`
4. `fluence_0p5_50ps_free_surface_ablation_stable_debug`
5. `fluence_0p5_100ps_free_surface_ablation_quickcheck`
6. `fluence_0p5_300ps_free_surface_ablation`

The first case must no longer show the `fix ttm/mod` non-periodic-boundary
error.  Do not run the 100 ps or 300 ps cases until the shorter free-surface
debug cases reach their final time without lost atoms or nonphysical max-KE /
max-vz spikes.  The final free-surface production case must reach 300 ps, not
30 ps.

## Expected Outputs

Each completed V8 case should be analyzed with `scripts/v80_postprocess_physical.py`.
The expected outputs are:

- `si_ablation_continuous_fields_early_v8.png`
- `si_ablation_continuous_fields_ejecta_v8.png`
- `overview_physical_v8.png`
- `surface_profile_evolution.png`
- `temperature_lattice_vs_time_v8.png`
- `pressure_vs_time_v8.png`
- `stability_diagnostics_v8.png`
- `depth_resolved_temperature_profile_v8.png`
- `electron_temperature_vs_time_v8.png` when `Te_out.*` files exist
- `ablation_metrics_physical_v8.csv`
- `surface_profiles_v8.csv`
- `heat_group_diagnostics_v8.csv`
- `stress_pressure_fields_v8.npz`
- `README_V8_physical_model.md`

If stable ejecta is still very small after 300 ps, the README should identify
whether the likely cause is low fluence, insufficient time, boundary effects,
source normalization, TTM parameters, or the potential's high-temperature
limitations.
