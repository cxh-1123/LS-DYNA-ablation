# V6 LAMMPS TTM-MD smoke summary

Pilot case `ep_5p0uj` (5.0 µJ, 30 ps pulse, fix ttm/mod, laser axis +x).

- profile: smoke (6144 atoms, 10 ps, dt=0.2 fs)
- wall time: ~24 min (12 OpenMP threads)
- lattice temp at 10 ps: ~4350 K
- snapshots: 0, 5, 10 ps
- full trajectories and figures stay in ignored `results/` / `figures/`

Regenerate inputs:

```powershell
python scripts\build_v6_lammps_ttm_input.py --smoke
powershell -ExecutionPolicy Bypass -File scripts\run_v6_selected.ps1 -OnlyCase ep_5p0uj -LammpsExe "D:\LAMMPS 64-bit 22Jul2025\bin\lmp.exe"
```
