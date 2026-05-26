# V3B plume / shock proxy calibration

V3B is a display and handoff proxy, not a calibrated plume CFD model.  Its
parameters should be treated as knobs until they are tied to experiment or
literature.

## What needs calibration

| Quantity | Current role | Calibration target |
| --- | --- | --- |
| `R_plume(t)` | plume front proxy | measured plume-front radius/height vs time |
| `R_shock(t)` | shock-front proxy | shadowgraph / schlieren shock front |
| density proxy | image contrast proxy | transmission or phase-contrast intensity |
| ejecta proxy | display-only particles | qualitative scatter / debris envelope |

## Minimal data to collect

- one early frame near `50-100 ps`;
- one frame near `0.5-1 ns`;
- one frame near `3-5 ns`;
- field of view scale in micrometres;
- pulse energy, spot radius, wavelength, and ambient condition.

## Suggested fitting order

1. Fix V2.6 driver geometry first: use the max crater-depth frame for `t0`,
   crater radius, and crater depth.
2. Fit plume-front growth with `v0_um_per_ns_power` and `b`.
3. Fit shock-front offset with `c0_um_per_ns`, `beta_um_per_sqrt_ns`, and
   `shock_ahead_min_um`.
4. Tune density contrast only after front positions look reasonable.
5. Tune ejecta particles last, and keep them labelled as display-only.

## Acceptance checks

- shock front should remain ahead of plume front;
- plume front should be monotonic after launch;
- density contrast should decay with time unless the experiment shows otherwise;
- ejecta particles must not be interpreted as conserved material mass;
- the final figures must state that V3B is a proxy model.

## Report wording

> The V3B plume, shock, and ejecta fields are calibrated proxy fields driven by
> the V2.6 threshold-equivalent crater metrics.  They are intended for
> visualization and structured-light handoff, not as a standalone multiphase
> ablation or gas-dynamics solution.
