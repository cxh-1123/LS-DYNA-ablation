# V7 diagnostics: fluence_0p5_100ps

This file is generated from the existing V61 LAMMPS trajectory.  No simulation
parameters were changed.

## Main diagnosis

- Current state: **early overheating / melting / expansion, with only weak initial ejecta**.
- Maximum surface atom temperature proxy p95: `1.15e+04 K`.
- Maximum vapor-like atom count by temperature proxy alone: `5591`.
- Maximum vapor-candidate count after density/height/velocity filters: `3291`.
- Maximum detached-candidate count: `7`.
- Maximum stable-ejecta count: `2`.
- Maximum robust p95 surface lift: `0.0493 nm`.
- Maximum raw max-z lift: `4.38 nm`.
- Raw lift outlier pollution: `yes`.
- Maximum low-density atoms with rho/rho0 < 0.3: `67`.
- Maximum pressure proxy p99: `44.3 GPa`.

## Why T > vapor can be large while stable ejecta is small

The per-atom temperature used here is a kinetic-energy proxy.  It can flag many
hot or strongly disordered atoms, but a real ejected atom must also be spatially
detached, low density, moving upward, and persistent across frames.  In this
trajectory, many atoms become hot, but most remain connected to the main slab.
That means the present 0.5 J/cm2, 100 ps run is not yet a strong plume case.

## Surface and crater metrics

`raw_lift_outlier_sensitive` uses the maximum z position and is kept only as a
diagnostic.  The main surface metric is `robust_lift`, computed from local
surface p95 values after excluding vapor/detached candidates and very low
density atoms.  Use `robust_lift`, `center_surface_z_p95`, `edge_surface_z_p95`,
`crater_depth_center`, and `rim_lift` for physical discussion.

## Pressure and stress

The pressure field is labelled `pressure_proxy`.  It is derived from per-atom
stress divided by an approximate initial atomic volume.  It is useful for
tracking stress-wave timing and sign, but should not be presented as a fully
calibrated hydrostatic pressure without additional validation.

## Is this enough for a paper-style main figure?

Not yet.  The continuous fields are useful diagnostics, but the current model is
too small laterally, only runs to 100 ps, and produces little stable ejecta.  V8
should update the physical model before treating the output as a main result.

## Likely V8 changes needed

- Increase total time to at least 300 ps.
- Use a larger lateral cell for pressure waves and plume formation.
- Re-check fluence normalization, reflected vs absorbed fluence, and optical
  penetration depth.
- Keep Si-specific TTM parameters explicit in the config.
- Add a controlled fluence scan rather than blindly increasing fluence.
- Treat Stillinger-Weber Si as a useful baseline but warn that high-temperature
  vapor/plume behavior requires validation.
