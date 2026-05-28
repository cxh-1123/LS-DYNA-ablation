V6 LAMMPS TTM-MD pilot -- ep_5p0uj
=====================================

Laser (aligned with V1.6/V1.7):
  Ep        = 5.0 uJ
  tau       = 30.0 ps
  w0        = 35.0 um
  A         = 0.5
  delta_abs = 5.0 nm

fix ttm/mod I_0 (metal units): 7.625733e-01 eV/(ps Ang^2)
Estimated peak volumetric scale (reference): 2.443e+21 W/m^3

Requirements:
  1. LAMMPS with EXTRA-FIX package (fix ttm/mod)
  2. potentials/Si.sw at repo root
  3. Run from this directory:
       lmp -in in_ep_5p0uj.lammps -log log_ep_5p0uj.lammps

Or:
  powershell -ExecutionPolicy Bypass -File scripts/run_v6_selected.ps1 -OnlyCase ep_5p0uj

Post-process:
  python scripts/plot_v6_lammps_ttm_snapshots.py --case ep_5p0uj
  python scripts/export_v6_density_for_v4.py --case ep_5p0uj

Outputs -> results/v6_lammps_ttm_md_pilot/ep_5p0uj/
