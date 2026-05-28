"""
build_v6_lammps_ttm_input.py
============================

Generate LAMMPS TTM-MD input decks for the V6 pilot (ep_5p0uj, 30 ps).

Uses fix ttm/mod (LAMMPS 22Jul2025 syntax) with laser along +x.
Writes under models/v6_lammps_ttm_md_pilot/<case>/:
  - in_<case>.lammps
  - <case>.ttm_mod
  - v6_run_manifest.json
  - README_run.txt

Does NOT run LAMMPS.  Does NOT modify V1.7 / V2.6 / V3B / V5A.

Run:
    python scripts\\build_v6_lammps_ttm_input.py
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import sys
import tomllib
from pathlib import Path

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


def load_toml(path: Path) -> dict:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def pulse_fluence_Wm2(Ep_J: float, w0_m: float, tau_s: float) -> float:
    """Peak surface intensity scale [W/m^2] for Gaussian temporal pulse."""
    return Ep_J / (math.sqrt(math.pi) * tau_s * math.pi * w0_m ** 2)


def peak_intensity_metal(Ep_J: float, w0_m: float, tau_s: float, A: float) -> float:
    """Peak absorbed intensity for fix ttm/mod I_0 [eV/(ps Ang^2)] in metal units."""
    i0_w_m2 = A * pulse_fluence_Wm2(Ep_J, w0_m, tau_s)
    return i0_w_m2 * 6.242e-14


def write_ttm_mod_file(
    path: Path,
    *,
    i0_metal: float,
    l_skin_ang: float,
    tau_ps: float,
    nx_grid: int,
) -> None:
    """
    fix ttm/mod init file (odd lines = comments).  Si coefficients from LAMMPS
    examples/ttm/Si.ttm_mod; laser I_0 / l_skin / tau from V1.6/V1.7 alignment.
    """
    lines = [
        "a_0, energy/(temperature*electron) units",
        "-0.00012899",
        "a_1, energy/(temperature^2*electron) units",
        "-0.0000000293276",
        "a_2, energy/(temperature^3*electron) units",
        "-0.0000229991",
        "a_3, energy/(temperature^4*electron) units",
        "-0.000000927036",
        "a_4, energy/(temperature^5*electron) units",
        "-0.0000011747",
        "C_0, energy/(temperature*electron) units",
        "0.000129",
        "A, 1/temperature units",
        "0.180501",
        "rho_e, electrons/volume units",
        "0.05",
        "D_e, length^2/time units",
        "20000",
        "gamma_p, mass/time units",
        "24.443",
        "gamma_s, mass/time units",
        "39.235",
        "v_0, length/time units",
        "79.76",
        "I_0, energy/(time*length^2) units",
        f"{i0_metal:.6e}",
        "lsurface, electron grid units (positive integer)",
        "0",
        "rsurface, electron grid units (positive integer)",
        str(nx_grid),
        "l_skin, length units",
        f"{l_skin_ang:.6g}",
        "tau, time units",
        f"{tau_ps:.6g}",
        "B, dimensionless",
        "0",
        "lambda, length units",
        "0",
        "n_ion, ions/volume units",
        "0.05",
        "surface_movement: 0 to disable tracking of surface motion, 1 to enable",
        "0",
        "T_e_min, temperature units",
        "300",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def apply_smoke_profile(cfg: dict) -> dict:
    smoke = cfg.get("profile", {}).get("smoke")
    if not smoke:
        return cfg
    out = copy.deepcopy(cfg)
    box = out["md_box"]
    box["nx_lat"] = int(smoke["nx_lat"])
    box["ny_lat"] = int(smoke["ny_lat"])
    box["nz_lat"] = int(smoke["nz_lat"])
    grid = out["ttm_grid"]
    grid["nx"] = int(smoke["ttm_nx"])
    grid["ny"] = int(smoke["ttm_ny"])
    grid["nz"] = int(smoke["ttm_nz"])
    sol = out["solver"]
    sol["t_end_ps"] = float(smoke["t_end_ps"])
    if "dt_fs" in smoke:
        sol["dt_fs"] = float(smoke["dt_fs"])
    sol["thermo_every"] = int(smoke["thermo_every"])
    sol["dump_every_ps"] = float(smoke["dump_every_ps"])
    sol["snapshot_times_ps"] = list(smoke["snapshot_times_ps"])
    out["meta"]["version"] = str(out["meta"].get("version", "V6.0")) + "+smoke"
    return out


def build_lammps_input(
    case: dict,
    cfg: dict,
    out_dir: Path,
    ttm_mod_name: str,
    i0_metal: float,
) -> Path:
    grid = cfg["ttm_grid"]
    sol = cfg["solver"]
    laser = cfg["laser"]
    pot = cfg["lammps"]["si_sw_potential"]
    pot_rel = Path(pot).as_posix()
    if not Path(pot).is_absolute():
        pot_rel = "../../../" + pot_rel
    dt_ps = float(sol["dt_fs"]) * 1e-3
    t_end = float(sol["t_end_ps"])
    nsteps = int(round(t_end / dt_ps))
    tau_ps = float(laser["pulse_width_ps"])
    snap = list(sol["snapshot_times_ps"])
    t0_k = float(cfg["material"]["T0_K"])

    in_path = out_dir / f"in_{case['name']}.lammps"

    snap_vars = "\n".join(
        f"variable tsnap{i} equal {t:g}" for i, t in enumerate(snap)
    )
    run_segments: list[str] = ["run             0"]
    prev_t = 0.0
    for t in snap:
        t_f = float(t)
        if t_f <= 0.0:
            run_segments.append(
                f'write_dump      all custom snapshot_{int(t_f)}ps.lammpstrj id type x y z vx vy vz'
            )
            continue
        nseg = int(round((t_f - prev_t) / dt_ps))
        if nseg > 0:
            run_segments.append(f"run             {nseg}")
        run_segments.append(
            f'write_dump      all custom snapshot_{int(t_f)}ps.lammpstrj id type x y z vx vy vz'
        )
        prev_t = t_f
    tail_steps = nsteps - int(round(prev_t / dt_ps))
    if tail_steps > 0:
        run_segments.append(f"run             {tail_steps}")
    run_block = "\n".join(run_segments)

    text = f"""# =============================================================================
# V6 LAMMPS TTM-MD pilot -- {case['name']}
# Generated by scripts/build_v6_lammps_ttm_input.py
# Ep = {case['Ep_uJ']} uJ, tau = {laser['pulse_width_ps']} ps, w0 = {laser['spot_radius_um']} um
# Laser axis: +x (fix ttm/mod).  NOT quantitative without TTM/Si.sw calibration.
# =============================================================================

units           metal
atom_style      atomic
boundary        {cfg['md_box']['boundary_x']} {cfg['md_box']['boundary_y']} {cfg['md_box']['boundary_z']}

variable        dt equal {dt_ps}
timestep        ${{dt}}

# --- geometry (diamond Si slab, depth along x) ---
lattice         diamond {cfg['material']['lattice_a_angstrom']}
variable        nxlat equal {cfg['md_box']['nx_lat']}
variable        nylat equal {cfg['md_box']['ny_lat']}
variable        nzlat equal {cfg['md_box']['nz_lat']}
region          sim block 0 ${{nxlat}} 0 ${{nylat}} 0 ${{nzlat}} units lattice
create_box      1 sim
create_atoms    1 box

mass            1 28.0855
pair_style      sw
pair_coeff      * * {pot_rel} Si

neighbor        1.0 bin
neigh_modify    every 1 delay 0 check yes

# --- TTM/mod (LAMMPS 22Jul2025+ syntax; laser along +x) ---
fix             twotemp all ttm/mod 1354684 {ttm_mod_name} {grid['nx']} {grid['ny']} {grid['nz']} set {t0_k}

# --- MD integrator (NVE + TTM coupling) ---
variable        xbot equal ${{nxlat}}-2
region          bottom block ${{xbot}} INF INF INF INF INF units lattice
group           frozen region bottom
fix             integ all nve
fix             freeze_bottom frozen setforce 0.0 0.0 0.0
velocity        all create {t0_k} 12345 rot yes dist gaussian
velocity        frozen set 0.0 0.0 0.0

# --- output ---
thermo_style    custom step time temp pe ke etotal f_twotemp[1] f_twotemp[2]
thermo          {int(sol['thermo_every'])}

dump            d1 all custom {max(int(round(float(sol['dump_every_ps'])/dt_ps)), 1)} traj.lammpstrj id type x y z vx vy vz
dump_modify     d1 sort id

{snap_vars}

# --- run (segmented for snapshot dumps) ---
{run_block}
write_data      final.data
"""
    in_path.write_text(text, encoding="utf-8")
    return in_path


def build_readme(
    out_dir: Path,
    case: dict,
    cfg: dict,
    i0_metal: float,
    peak_Wm3: float,
) -> None:
    txt = f"""V6 LAMMPS TTM-MD pilot -- {case['name']}
=====================================

Laser (aligned with V1.6/V1.7):
  Ep        = {case['Ep_uJ']} uJ
  tau       = {cfg['laser']['pulse_width_ps']} ps
  w0        = {cfg['laser']['spot_radius_um']} um
  A         = {cfg['laser']['absorption_A']}
  delta_abs = {cfg['laser']['absorption_depth_nm']} nm

fix ttm/mod I_0 (metal units): {i0_metal:.6e} eV/(ps Ang^2)
Estimated peak volumetric scale (reference): {peak_Wm3:.3e} W/m^3

Requirements:
  1. LAMMPS with EXTRA-FIX package (fix ttm/mod)
  2. potentials/Si.sw at repo root
  3. Run from this directory:
       lmp -in in_{case['name']}.lammps -log log_{case['name']}.lammps

Or:
  powershell -ExecutionPolicy Bypass -File scripts/run_v6_selected.ps1 -OnlyCase {case['name']}

Post-process:
  python scripts/plot_v6_lammps_ttm_snapshots.py --case {case['name']}
  python scripts/export_v6_density_for_v4.py --case {case['name']}

Outputs -> results/v6_lammps_ttm_md_pilot/{case['name']}/
"""
    (out_dir / "README_run.txt").write_text(txt, encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--toml", type=Path, default=root / "config/v6_lammps_ttm_md_pilot.toml")
    ap.add_argument("--out-root", type=Path, default=root / "models/v6_lammps_ttm_md_pilot")
    ap.add_argument("--smoke", action="store_true", help="Use profile.smoke (small box, 50 ps)")
    args = ap.parse_args()

    cfg = load_toml(args.toml)
    if args.smoke:
        cfg = apply_smoke_profile(cfg)
        print("[INFO] smoke profile: smaller box + t_end=50 ps")
    args.out_root.mkdir(parents=True, exist_ok=True)
    (root / "potentials").mkdir(parents=True, exist_ok=True)

    laser = cfg["laser"]
    w0_m = float(laser["spot_radius_um"]) * 1e-6
    tau_s = float(laser["pulse_width_ps"]) * 1e-12
    delta_m = float(laser["absorption_depth_nm"]) * 1e-9
    A = float(laser["absorption_A"])
    a_lat = float(cfg["material"]["lattice_a_angstrom"])
    lx = cfg["md_box"]["nx_lat"] * a_lat
    ly = cfg["md_box"]["ny_lat"] * a_lat
    lz = cfg["md_box"]["nz_lat"] * a_lat
    nx_grid = int(cfg["ttm_grid"]["nx"])
    l_skin_ang = delta_m * 1e10
    tau_ps = float(laser["pulse_width_ps"])

    registry: list[dict[str, str]] = []

    for case in cfg["cases"]:
        name = str(case["name"])
        Ep_J = float(case["Ep_uJ"]) * 1e-6
        peak_flux = pulse_fluence_Wm2(Ep_J, w0_m, tau_s)
        peak_Wm3 = A * peak_flux / max(delta_m, 1e-12)
        i0_metal = peak_intensity_metal(Ep_J, w0_m, tau_s, A)

        out_dir = args.out_root / name
        out_dir.mkdir(parents=True, exist_ok=True)

        ttm_mod_name = f"{name}.ttm_mod"
        write_ttm_mod_file(
            out_dir / ttm_mod_name,
            i0_metal=i0_metal,
            l_skin_ang=l_skin_ang,
            tau_ps=tau_ps,
            nx_grid=nx_grid,
        )

        in_path = build_lammps_input(case, cfg, out_dir, ttm_mod_name, i0_metal)
        build_readme(out_dir, case, cfg, i0_metal, peak_Wm3)

        manifest = {
            "case": name,
            "profile": "smoke" if args.smoke else "full",
            "Ep_uJ": case["Ep_uJ"],
            "pulse_width_ps": laser["pulse_width_ps"],
            "spot_radius_um": laser["spot_radius_um"],
            "I0_metal_eV_ps_A2": i0_metal,
            "peak_Wm3_estimate": peak_Wm3,
            "slab_angstrom": {"lx": lx, "ly": ly, "lz": lz},
            "laser_axis": "+x",
            "ttm_style": "ttm/mod",
            "snapshot_times_ps": cfg["solver"]["snapshot_times_ps"],
            "input": str(in_path.relative_to(root)).replace("\\", "/"),
        }
        (out_dir / "v6_run_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        registry.append({
            "name": name,
            "Ep_uJ": str(case["Ep_uJ"]),
            "input_lammps": f"models/v6_lammps_ttm_md_pilot/{name}/in_{name}.lammps",
            "run_lammps": str(case.get("run_lammps", True)).lower(),
            "t_end_ps": str(cfg["solver"]["t_end_ps"]),
            "notes": str(case.get("notes", "")),
        })
        print(f"[OK] {name}: in={in_path.name}, ttm_mod={ttm_mod_name}, I0~{i0_metal:.3e}")

    reg_path = args.out_root / "v6_case_registry.csv"
    with reg_path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(registry[0].keys()))
        wr.writeheader()
        wr.writerows(registry)
    print(f"[OK] registry -> {reg_path}")
    print("\nNext: copy Si.sw to potentials/Si.sw, then run scripts/run_v6_selected.ps1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
