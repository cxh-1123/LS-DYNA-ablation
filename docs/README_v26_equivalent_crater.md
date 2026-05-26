# V2.6 threshold-equivalent crater definition

V2.6 does not remove LS-DYNA elements and does not simulate real material
ejection. It reads the V1.7 temperature history and applies temperature
thresholds:

- melt-equivalent zone: `T >= T_melt = 1687 K`
- vapor / ablation-equivalent zone: `T >= T_vap = 3538 K`
- equivalent crater depth: vapor-threshold depth on the symmetry axis
- equivalent crater radius: vapor-threshold radius on the top surface

This is useful because it converts a temperature field into simple, comparable
numbers. It is not a physical crater with mass conservation, recoil pressure,
phase explosion, stress wave fracture, or SPH/ALE ejecta.

Recommended wording for reports:

> The V2.6 crater metrics are temperature-threshold-equivalent quantities
> derived from V1.7 temperature fields. The reported crater depth and radius
> indicate where the lattice temperature exceeds the vapor threshold; they do
> not represent explicit LS-DYNA material removal or a fully resolved
> multiphase ablation process.

Recommended validation checks:

- peak temperature should stay at the axis/top surface;
- melt/vapor radius should shrink as the threshold increases;
- crater depth and radius should be compared across energy and mesh resolution;
- the same `T_melt` and `T_vap` values must be used in all V26 plots and CSVs.
