# V1.7 grid convergence

This workflow prepares a focused surface-mesh sensitivity check for the 30 ps
local thermal model.  It compares surface `dz = 5 / 10 / 20 nm` while keeping
the same laser, material, time window, and radial grid.

Default comparison cases:

- `ep_1p5uj`: melt-onset bracket
- `ep_2p5uj`: vapor-onset bracket
- `ep_5p0uj`: high-energy ablation candidate

Build inputs only:

```powershell
python scripts\build_v17_grid_convergence.py --dry-run
python scripts\build_v17_grid_convergence.py
```

Suggested acceptance rule after running LS-DYNA:

- compare `T_peak_K`, max melt depth/radius, and max vapor depth/radius;
- if the `10 nm` result differs from `5 nm` by less than about `5-10%`,
  the current V17 mesh is acceptable;
- if `20 nm` differs strongly while `10 nm` and `5 nm` agree, keep `10 nm`;
- if `10 nm` and `5 nm` still differ strongly, refine the surface mesh or
  shorten the thermal timestep before trusting threshold numbers.

This is still a single-temperature LS-DYNA thermal comparison, not a full TTM
or explicit material-removal validation.
