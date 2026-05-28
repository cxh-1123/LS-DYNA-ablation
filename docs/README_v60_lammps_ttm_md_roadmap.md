# V6.0 LAMMPS TTM-MD Roadmap for Explicit Ejecta/Plume Evolution

## Goal

Build a LAMMPS-based TTM-MD workflow that can produce paper-style laser ablation figures:

- atomistic ejecta/plume evolution over time;
- lattice temperature field and molten/solid interface;
- pressure/stress wave evolution inside silicon;
- surface deformation and ablation crater morphology;
- unified figure style compatible with the existing V45/V50 panels.

This replaces the V50 "explicit ejecta marker layer" proxy with a physically resolved atomistic ejecta process.  V50 remains useful as the visual target and plotting prototype.

## Physical Model Choice

Use TTM-MD:

- atoms are evolved by molecular dynamics;
- electrons are represented by a continuum temperature grid;
- laser energy is deposited into the electronic subsystem;
- electron energy diffuses and transfers to atoms through electron-phonon coupling;
- ablation/ejecta are not manually prescribed: atoms leave the surface when the MD dynamics provide enough energy and momentum.

LAMMPS supports this through `fix ttm`, `fix ttm/grid`, `fix ttm/mod`, and `fix ttm/thermal`.  The first robust path should use `fix ttm` or `fix ttm/mod`; `fix ttm/thermal` is promising but should only be used after confirming the local LAMMPS build includes a sufficiently recent version.

## Stage V6A: Small, Verifiable TTM-MD Silicon Slab

Purpose: prove that LAMMPS, the Si potential, the TTM fix, and the dump pipeline work.

Model:

- material: crystalline Si;
- geometry: nanoscale slab with vacuum above the surface;
- boundary: periodic in lateral directions, fixed or thermostatted bottom layer, free top surface;
- potential: start with Stillinger-Weber Si (`pair_style sw`, `si.sw`) for robust melting/ablation behavior; keep Tersoff as a sensitivity check;
- timestep: about 0.1-0.25 fs for high-temperature ablation runs;
- simulation window: 0-100 ps first, then extend if stable.

Outputs:

- atom dump every 0.5-2 ps;
- per-atom velocity, kinetic energy, potential energy;
- local atomic temperature by depth bins;
- electronic temperature grid from TTM output;
- total energy audit.

Pass criteria:

- stable equilibration at 300 K;
- laser energy produces reasonable surface heating;
- energy balance does not drift catastrophically;
- no artificial atom loss before irradiation.

## Stage V6B: Laser Source Calibration

Purpose: connect the LAMMPS energy input to the LS-DYNA/V17/V26 threshold work.

Inputs to calibrate:

- pulse width: 30 ps baseline;
- absorbed fluence/energy: map 5 uJ and spot size to absorbed fluence;
- optical penetration depth or source depth;
- electron heat capacity, electron thermal conductivity, and electron-phonon coupling.

Calibration targets:

- onset of melting near the previous threshold range;
- onset of material removal near the V26/V41 5 uJ case;
- lattice temperature peak and cooling trend comparable to the earlier thermal model, but not forced to match exactly.

Important note:

The LS-DYNA single-temperature model sends energy directly to the lattice.  TTM-MD sends energy first to electrons and then to atoms, so early-time lattice temperature should lag behind electron temperature.

## Stage V6C: Explicit Atomistic Ejecta

Purpose: generate the plume/ejecta from actual atom motion.

Definitions:

- substrate atoms: atoms below the moving surface;
- ejecta atoms: atoms above the original surface with upward velocity and separation from the connected solid;
- molten atoms: atoms with local order loss or local lattice temperature above melting threshold;
- vapor/plume atoms: atoms above the surface with high kinetic energy and low local coordination.

Analysis outputs:

- plume snapshots at 5, 15, 25, 40, 55, 70, 85, 95 ps;
- atom color by lattice temperature or kinetic temperature;
- atom size/opacity separated by solid, liquid, vapor/ejecta classes;
- plume height and front velocity versus time;
- number of emitted atoms versus time;
- ablation depth/crater depth versus time.

## Stage V6D: Temperature, Pressure, and Phase Panels

Purpose: produce figures comparable to the paper examples.

Temperature:

- bin atoms by x-z or r-z section;
- compute lattice temperature from thermal velocity after subtracting local center-of-mass velocity;
- overlay electronic temperature from TTM grid if available.

Pressure/stress:

- compute per-atom stress in LAMMPS;
- bin stress into spatial cells;
- pressure = negative one third of the trace of the local stress tensor;
- plot pressure wave inside the silicon with the same palette and layout as V45/V50.

Phase:

- use local order/coordination plus temperature threshold to classify solid, liquid, vapor/ejecta;
- show liquid-crystal interface as a contour, not only colored atoms.

## Stage V6E: Unified Figure Style

All V6 figures should use one consistent visual grammar:

- same time labels: ps for <1 ns, ns above that;
- same axes: signed lateral coordinate centered on laser spot, depth/height vertical;
- same color maps:
  - temperature: purple/blue to yellow/red, log or clipped linear as needed;
  - pressure: diverging map centered at 0 GPa;
  - displacement/height: sequential map;
  - ejecta atoms: bright warm colors over a light background;
- same panel ordering:
  - row 1: atomistic plume/ejecta snapshots;
  - row 2: lattice temperature or phase map;
  - row 3: pressure/stress wave;
  - side panel: peak temperature, emitted atom count, plume front, crater depth.

## Stage V6F: Validation and Sensitivity

Minimum checks before treating results as physically meaningful:

- energy conservation audit with atomic + electronic energy;
- timestep convergence: 0.25 fs vs 0.1 fs;
- lateral size convergence: plume and crater not dominated by periodic images;
- depth convergence: bottom boundary does not reflect pressure waves too early;
- TTM grid convergence in depth and lateral direction;
- potential sensitivity: SW vs Tersoff if feasible;
- fluence sweep around threshold: below ablation, near threshold, above threshold.

## Practical Build Order

1. Add `config/v60_lammps_ttm_md_pilot.toml`.
2. Add `scripts/build_v60_lammps_ttm_md_pilot.py` to generate:
   - LAMMPS data file;
   - LAMMPS input script;
   - TTM initial electron temperature file;
   - run registry.
3. Add `scripts/run_v60_selected.ps1`.
4. Run a tiny dry test: no laser, 2 ps.
5. Run low-fluence heating: no ablation expected.
6. Run 5 uJ-equivalent pilot: look for ejecta atoms.
7. Add `scripts/analyze_v60_lammps_dump.py`.
8. Add `scripts/plot_v60_ttm_md_evolution.py` using the V45/V50 style.
9. Push only lightweight outputs:
   - CSV summaries;
   - final PNG panels;
   - input decks/configs/scripts;
   - not huge dump trajectories.

## Immediate Next Step

Start with V6A, not the full paper figure.  The first runnable target is:

> a small Si slab TTM-MD pilot that outputs atom snapshots, electronic temperature grid, lattice temperature by depth, and an energy audit.

Once V6A is stable, V6B/V6C can make the laser source strong enough to produce atomistic ejection.

## References to Check and Cite

- LAMMPS `fix ttm`, `fix ttm/grid`, `fix ttm/mod`, and `fix ttm/thermal` documentation: https://docs.lammps.org/fix_ttm.html
- LAMMPS `pair_style sw` documentation for Stillinger-Weber Si: https://docs.lammps.org/pair_sw.html
- LAMMPS `pair_style tersoff` documentation for Si Tersoff variants: https://docs.lammps.org/pair_tersoff.html
- LAMMPS dump/visualization documentation: https://docs.lammps.org/dump.html
- Duffy and Rutherford, J. Phys.: Condens. Matter 19, 016207 (2007): https://doi.org/10.1088/0953-8984/19/1/016207
