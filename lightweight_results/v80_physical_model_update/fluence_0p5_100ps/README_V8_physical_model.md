# README V8 physical model diagnostics: fluence_0p5_100ps

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

- Initial mean bin density check: `1`; this should be close
  to 1 for the density normalization to be usable.
- Maximum low-density area fraction (`rho/rho0 < 0.7`): `0.258`.
- Maximum void area fraction (`rho/rho0 < 0.3`): `0.0537`.
- Maximum raw max-z lift: `4.41 nm`.
- Maximum robust p95 lift: `0.409 nm`.
- Raw-lift outlier warning: `YES`.
- Maximum detached candidates: `2`.
- Maximum stable ejecta: `2`.
- Maximum pressure p99 after bin conversion: `15.3 GPa`.

## Interpretation

The existing `0.5 J/cm2, 100 ps` run remains best interpreted as early
overheating, lattice disordering and weak atom escape.  It is not yet a stable
plume result because the largest connected component still contains nearly all
atoms and stable ejecta remains very small.

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
