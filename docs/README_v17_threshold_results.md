# V1.7 30 ps threshold result interpretation

This note interprets the lightweight summary in
`lightweight_results/v17_30ps_local/v17_case_summary.csv`.

## Key thresholds

Material reference temperatures:

- melt temperature: `T_melt = 1687 K`
- vapor temperature: `T_vap = 3538 K`

Observed V1.7 single-temperature results:

| Case | Pulse energy (uJ) | Peak temperature (K) | Regime |
| --- | ---: | ---: | --- |
| `ep_1p0uj` | 1.00 | 1361.60 | no melt |
| `ep_1p5uj` | 1.50 | 1891.87 | melt only |
| `ep_3p0uj` | 3.00 | 3483.53 | melt only |
| `ep_3p05uj` | 3.05 | 3535.87 | melt only, just below vapor |
| `ep_3p1uj` | 3.10 | 3590.86 | vapor/ablation candidate |
| `ep_3p2uj` | 3.20 | 3696.44 | vapor/ablation candidate |

Plain-language conclusion:

- The melt threshold is between `1.0` and `1.5 uJ`.
- Linear interpolation gives an estimated melt onset near `1.31 uJ`.
- The vapor threshold is between `3.05` and `3.10 uJ`.
- Linear interpolation gives an estimated vapor onset near `3.052 uJ`.

## Interpretation

The new `3.05 / 3.1 / 3.2 uJ` cases make the vapor threshold much clearer.
The `3.05 uJ` case peaks at `3535.87 K`, only about `2.13 K` below the
nominal vapor temperature.  The `3.1 uJ` case reaches `3590.86 K`, so it is
already above the vapor-temperature criterion.

For reports, the safest wording is:

> In the V1.7 30 ps single-temperature LS-DYNA approximation, the lattice
> vapor-threshold energy is bracketed between 3.05 and 3.10 uJ, with a simple
> linear interpolation estimate of about 3.052 uJ.  This should be interpreted
> as a temperature-threshold estimate, not as explicit material removal.

## Remaining useful checks

The current lightweight result package still has missing runs for
`ep_1p2uj`, `ep_1p8uj`, and `ep_2p2uj`.  These are not required to identify the
vapor threshold, but `ep_1p2uj` would help confirm the melt-threshold bracket.

For a tighter melt threshold, useful next cases would be around `1.30` and
`1.35 uJ`.
