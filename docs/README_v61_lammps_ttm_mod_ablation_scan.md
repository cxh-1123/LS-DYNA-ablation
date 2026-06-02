# V6B LAMMPS `ttm/mod` Ablation Scout

V6B moves the project in the physically correct direction for ultrafast silicon
ablation:

1. 0-1 ps: electrons absorb laser energy first.
2. 1-10 ps: electron-phonon coupling heats the lattice.
3. 10-30 ps: the 30 ps pulse continues depositing energy; melting and stress
   waves may begin.
4. 30-60 ps: electron temperature falls while lattice heating, expansion, and
   pressure-wave propagation continue.
5. 60-100 ps: ejecta may appear if the fluence is above threshold.

The current V6A result was useful but weak: electron temperature existed, the
lattice only warmed modestly, surface lift was about 0.04 nm, and ejecta count
was zero.  V6B therefore switches from an initial hot-electron condition to the
LAMMPS `fix ttm/mod` laser-source workflow.

## Cases

The first stable scout cases are:

- `fluence_0p2_100ps`
- `fluence_0p5_100ps`
- `fluence_1p0_100ps`

These are not final calibrated results.  They are threshold-finding cases to
identify whether the current Si/SW/TTM parameter set can produce melting and
ejecta within 100 ps.

## Run Order

Build:

```powershell
python scripts\build_v61_lammps_ttm_mod_ablation_scan.py
```

Dry run:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v61_selected.ps1 -DryRun
```

Run the middle case first:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v61_selected.ps1 `
  -OnlyCase fluence_0p5_100ps `
  -LammpsExe "D:\LAMMPS 64-bit 22Jul2025\bin\lmp.exe" `
  -PotentialFile "D:\cxh-daima\LS-DYNA-ablation\potentials\Si.sw"
```

If it stays too weak, run `fluence_1p0_100ps`.  If it becomes unstable, run
`fluence_0p2_100ps`.

## Interpretation

The goal is not just to see a plume.  A physically credible result should show:

- high electron temperature before lattice temperature response;
- lattice temperature crossing the Si melting point near the surface;
- stress/pressure-wave formation after rapid heating;
- surface expansion followed by separated high-velocity atoms or clusters;
- ejecta timing in the tens to hundreds of ps range, not immediately at time
  zero.

## Reference Basis

- LAMMPS example `Examples/ttm/in.ttm.mod.lmp` and `Examples/ttm/Si.ttm_mod`
  from the local LAMMPS 22Jul2025 install.
- LAMMPS `fix ttm/mod` citation printed by LAMMPS:
  Pisarev et al., J. Phys.: Condens. Matter 26, 475401 (2014);
  Norman et al., Contrib. Plasma Phys. 53, 129-139 (2013).
