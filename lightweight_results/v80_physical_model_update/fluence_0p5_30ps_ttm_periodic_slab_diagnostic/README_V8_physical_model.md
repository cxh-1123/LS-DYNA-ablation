# README V8 physical model diagnostics: fluence_0p5_30ps_ttm_periodic_slab_diagnostic

## Case scope

This is a TTM-MD periodic slab diagnostic, not strict free-surface ablation.  Use it only to verify electron-lattice coupling and temperature evolution.

This V8 diagnostic pass is generated from the existing trajectory.  It fixes
the V7 postprocessing definitions before changing the simulation itself.

## V7 issues addressed

- Density now uses `bin_density = current x-z bin count / initial x-z bin count`
  for continuous fields.  Neighbor density is kept only for connectivity and
  ejecta classification.
- Atom temperature now reports both raw kinetic temperature and
  local-drift-subtracted thermal temperature.  Main figures use thermal
  temperature.
- Pressure and von Mises stress use bin-summed LAMMPS `stress/atom` divided by
  bin volume.  LAMMPS documents `compute stress/atom` as pressure*volume, so
  division by a meaningful volume is required before interpreting it as stress.
- Surface lift uses robust p95 profiles after excluding detached and
  low-neighbor-density atoms.

## Current trajectory diagnosis

- Heated surface region definition: atoms from `3 nm`
  below the initial surface to `0.5 nm` above it,
  after excluding detached/low-density candidates.  Atom count range in this
  region: `16354` to `17498`.
- Continuous temperature/density bins require at least
  `2` atoms; stress/pressure bins require at least
  `3` atoms.
- Initial mean bin density check: `1`; this should be close
  to 1 for the density normalization to be usable.
- Maximum low-density area fraction (`rho/rho0 < 0.7`): `0.273`.
- Maximum void area fraction (`rho/rho0 < 0.3`): `0.0856`.
- Maximum raw max-z lift: `0.904 nm`.
- Maximum robust p95 lift: `0.374 nm`.
- Raw-lift outlier warning: `limited`.
- Maximum detached candidates: `0`.
- Maximum stable ejecta: `0`.
- Maximum pressure p99 after bin conversion: `5.6 GPa`.
- Maximum global lattice/atom temperature from thermo: `2.14e+03 K`.
- Maximum heated-surface thermal temperature p95: `4.31e+03 K`.
- Maximum heated-surface thermal temperature max is retained only as a hotspot
  reference.  Do not use the max curve as the primary physical conclusion.
- Maximum molten atom fraction by thermal temperature: `0.474`.
- Maximum vapor-reference atom fraction by thermal temperature: `0.161`.

## TTM heat/source distribution check

`v80_fluence_0p5_30ps_ttm_periodic_slab_diagnostic.ttm_mod` reports `I_0 = 0.0127043`, `l_skin = 100 A` (10 nm), `tau = 30 ps`, and `D_e = 20000`.  The source is not intended as uniform slab heating, but the current penetration depth and electron diffusion parameters can still make the 30 ps diagnostic look too deep or too spatially uniform.  Treat this as a parameter-validation warning before any paper-style Si interpretation.

## Interpretation

Interpret this case according to its declared scope above.  If stable ejecta is
small while raw max-z lift is much larger than robust p95 lift, the apparent
surface motion is dominated by a few high-z atoms rather than a coherent crater
or plume.

For the periodic-slab diagnostic case, do not use these figures to claim a
real crater, stable ejecta, plume, or final morphology.  Its only valid role is
checking whether the TTM-MD coupling, electron temperature output, early lattice
heating, and melt-threshold crossing behave plausibly.

Vapor, detached and stable-ejecta labels in the overview are algorithmic
candidates.  For a periodic slab diagnostic, they are not physical plume or
ejecta conclusions.

## Pressure/stress status

Pressure/stress are now converted using LAMMPS metal-unit virial stress divided
by x-z-bin volume and expressed in GPa.  They are more physical than V7's single
atom-volume proxy, but still require caution because the bin volume is a
coarse-grained control volume and low-count bins are masked.

## What V8 simulation must change next

- The laser source must distinguish incident and absorbed fluence using Si
  reflectivity and optical penetration depth.
- A 30 ps FWHM pulse must be used; 100 ps is an observation window, not a pulse
  duration.
- Lateral size must increase beyond the debug ~4 nm width.
- Run time must increase to 300 ps before judging stable ejecta and plume
  formation.
- Fluence scan should start with 0.5, 0.8 and 1.0 J/cm2, then expand only if the
  simulations remain stable.
