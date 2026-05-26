# V1.7 single-temperature limitation

V1.7 is a single-temperature LS-DYNA thermal model.  At 30 ps this is a useful
engineering approximation, but it is not a full picosecond laser-ablation
model.

Plain-language interpretation:

- the laser energy is deposited directly into the LS-DYNA thermal field;
- real silicon first stores a large part of the energy in excited electrons;
- the lattice temperature follows after electron-phonon coupling over roughly
  tens of picoseconds;
- therefore V1.7 should be compared with the lattice branch `Tl` from the V1.6
  TTM prototype, not with electron temperature `Te`.

Use V1.7 for:

- checking whether the refined 30 ps mesh behaves numerically;
- estimating lattice-scale melt/vapor threshold brackets;
- producing temperature fields for threshold-equivalent V2.6 post-processing.

Do not use V1.7 alone for:

- quantitative ultrafast carrier dynamics;
- calibrated electron-lattice nonequilibrium;
- phase explosion, stress confinement, or material ejection;
- real ablation mass conservation.

Recommended wording for reports:

> The V1.7 model is a local refined 2D axisymmetric single-temperature
> LS-DYNA thermal approximation.  It resolves the near-surface thermal gradient
> for a 30 ps pulse, but it does not explicitly solve two-temperature
> electron-lattice dynamics.  Thresholds derived from V1.7 should therefore be
> interpreted as lattice-temperature, single-temperature estimates and checked
> against the V1.6 TTM prototype trends.
