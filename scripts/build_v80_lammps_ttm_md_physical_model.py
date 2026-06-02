"""
Build V8 Si laser-ablation LAMMPS inputs.

V8 deliberately separates two different physical meanings:

1. ttm_periodic_slab_diagnostic
   - boundary p p p
   - uses fix ttm/mod
   - short 30 ps diagnostic only

2. free_surface_ablation
   - boundary p p fs
   - does not use fix ttm/mod
   - uses an equivalent Beer-Lambert lattice heat source
   - main morphology/ejecta/plume visualization case
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore


if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


J_CM2_TO_EV_A2 = 624.1509074460763


def load_toml(path: Path) -> dict:
    with path.open("rb") as f:
        return tomllib.load(f)


def steps(ps: float, dt: float) -> int:
    return max(1, int(round(ps / dt)))


def run_value(cfg: dict, case: dict, key: str):
    return case.get(key, cfg["run_control"].get(key))


def case_timestep_ps(cfg: dict, case: dict) -> float:
    return float(run_value(cfg, case, "timestep_ps"))


def gaussian_sigma_ps(fwhm_ps: float) -> float:
    return fwhm_ps / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def absorbed_fluence(cfg: dict, case: dict) -> float:
    return float(case["incident_fluence_J_cm2"]) * (1.0 - float(cfg["laser"]["reflectivity_R"]))


def gaussian_peak_intensity(absorbed_j_cm2: float, fwhm_ps: float, scale: float) -> float:
    sigma = gaussian_sigma_ps(fwhm_ps)
    absorbed_eV_A2 = absorbed_j_cm2 * J_CM2_TO_EV_A2
    return scale * absorbed_eV_A2 / (sigma * math.sqrt(2.0 * math.pi))


def ttm_text(cfg: dict, case: dict) -> str:
    p = cfg["ttm_parameters_si"]
    laser = cfg["laser"]
    absorbed = absorbed_fluence(cfg, case)
    I0 = gaussian_peak_intensity(absorbed, float(laser["pulse_fwhm_ps"]), float(laser["intensity_scale"]))
    values = [
        ("a_0, energy/(temperature*electron) units", p["a_0"]),
        ("a_1, energy/(temperature^2*electron) units", p["a_1"]),
        ("a_2, energy/(temperature^3*electron) units", p["a_2"]),
        ("a_3, energy/(temperature^4*electron) units", p["a_3"]),
        ("a_4, energy/(temperature^5*electron) units", p["a_4"]),
        ("C_0, energy/(temperature*electron) units", p["C_0"]),
        ("A, 1/temperature units", p["A"]),
        ("rho_e, electrons/volume units", p["rho_e"]),
        ("D_e, length^2/time units", p["D_e"]),
        ("gamma_p, mass/time units", p["gamma_p"]),
        ("gamma_s, mass/time units", p["gamma_s"]),
        ("v_0, length/time units", p["v_0"]),
        ("I_0, energy/(time*length^2) units", I0),
        ("lsurface, electron grid units (positive integer)", 0),
        ("rsurface, electron grid units (positive integer)", 1),
        ("l_skin, length units", float(laser["optical_penetration_depth_nm"]) * 10.0),
        ("tau, time units", float(laser["pulse_fwhm_ps"])),
        ("B, dimensionless", p["B"]),
        ("lambda, length units", p["lambda"]),
        ("n_ion, ions/volume units", p["n_ion"]),
        ("surface_movement: 0 to disable tracking of surface motion, 1 to enable", p["surface_movement"]),
        ("T_e_min, temperature units", p["T_e_min"]),
    ]
    lines: list[str] = []
    for label, value in values:
        lines.append(str(label))
        lines.append(f"{float(value):.12g}" if isinstance(value, float) else str(value))
    return "\n".join(lines) + "\n"


def common_preamble(cfg: dict, case: dict, boundary: str) -> tuple[str, dict[str, float]]:
    size = cfg["sizes"][case["size"]]
    mat = cfg["material"]
    bnd = cfg["boundary"]
    run = cfg["run_control"]
    dt = case_timestep_ps(cfg, case)
    nx, ny, nz = int(size["nx_cells"]), int(size["ny_cells"]), int(size["nz_cells"])
    nz_box = nz + int(size["vacuum_cells_z"])
    alat = float(size["lattice_constant_A"])
    top_A = nz * alat
    area_A2 = (nx * alat) * (ny * alat)
    mobile_lo = float(bnd["bottom_fixed_thickness_A"])
    bath_hi = mobile_lo + float(bnd["bottom_thermostat_thickness_A"])
    eq_steps = steps(float(run["equilibration_ps"]), dt)
    thermo_steps = int(run_value(cfg, case, "thermo_interval_steps"))
    text = f"""units           metal
dimension       3
boundary        {boundary}
atom_style      atomic

variable        alat equal {alat:.6f}
lattice         diamond ${{alat}}
region          box block 0 {nx} 0 {ny} 0 {nz_box} units lattice
create_box      1 box
region          si block 0 {nx} 0 {ny} 0 {nz} units lattice
create_atoms    1 region si

mass            1 {float(mat['mass_amu']):.6f}
pair_style      {mat['pair_style']}
pair_coeff      * * {mat['pair_coeff_file']} Si

neighbor        2.0 bin
neigh_modify    every 5 delay 0 check yes

region          bottom block INF INF INF INF INF {mobile_lo:.6f} units box
region          bath block INF INF INF INF {mobile_lo:.6f} {bath_hi:.6f} units box
group           bottom region bottom
group           bath region bath
group           mobile subtract all bottom

compute         peatom all pe/atom
compute         keatom all ke/atom
compute         satom all stress/atom NULL

timestep        {dt:.8f}
velocity        mobile create 300.0 4928459 mom yes rot yes dist gaussian
velocity        bottom set 0.0 0.0 0.0
fix             hold bottom setforce 0.0 0.0 0.0

min_style       cg
minimize        1.0e-8 1.0e-10 500 5000

fix             nve mobile nve
fix             bathfix bath langevin 300.0 300.0 {float(bnd['bottom_langevin_damping_ps']):.6f} 827364 zero yes
thermo          {thermo_steps}
thermo_style    custom step time temp pe ke etotal press
run             {eq_steps}
unfix           bathfix
reset_timestep  0
"""
    meta = {
        "dt": dt,
        "top_A": top_A,
        "area_A2": area_A2,
        "eq_steps": eq_steps,
        "thermo_steps": thermo_steps,
    }
    return text, meta


def dump_and_run_lines(case: dict, cfg: dict, dump_prefix: str, prod_steps: int) -> str:
    run = cfg["run_control"]
    dt = case_timestep_ps(cfg, case)
    early_steps = min(steps(float(run["early_window_end_ps"]), dt), prod_steps)
    late_steps = max(0, prod_steps - early_steps)
    dump_early = steps(float(run_value(cfg, case, "dump_interval_early_ps")), dt)
    dump_late = steps(float(run_value(cfg, case, "dump_interval_late_ps")), dt)
    restart_steps = steps(float(run_value(cfg, case, "restart_interval_ps") or 10.0), dt)
    name = case["name"]
    return f"""
variable        vz_atom atom vz
variable        z_atom atom z
compute         max_ke all reduce max c_keatom
compute         max_vz all reduce max v_vz_atom
compute         zmax all reduce max v_z_atom
thermo          {int(run_value(cfg, case, "thermo_interval_steps"))}
thermo_style    custom step time temp pe ke etotal press c_max_ke c_max_vz c_zmax
restart         {restart_steps} restart.v80_{name}.1 restart.v80_{name}.2
dump            atomdump_early all custom {dump_early} dump.v80_{name}.early.lammpstrj id type x y z vx vy vz c_keatom c_peatom c_satom[1] c_satom[2] c_satom[3] c_satom[4] c_satom[5] c_satom[6]
dump_modify     atomdump_early sort id
run             {early_steps}
undump          atomdump_early
dump            atomdump_late all custom {dump_late} dump.v80_{name}.late.lammpstrj id type x y z vx vy vz c_keatom c_peatom c_satom[1] c_satom[2] c_satom[3] c_satom[4] c_satom[5] c_satom[6]
dump_modify     atomdump_late sort id
run             {late_steps}
undump          atomdump_late

write_data      final.v80_{name}.data

uncompute       max_ke
uncompute       max_vz
uncompute       zmax
unfix           nve
unfix           hold
"""


def ttm_diagnostic_input(cfg: dict, case: dict) -> str:
    grid = cfg["ttm_grid"]
    name = case["name"]
    dt = case_timestep_ps(cfg, case)
    prod_steps = steps(float(case["production_ps"]), dt)
    preamble, _ = common_preamble(cfg, case, boundary="p p p")
    dump_early = steps(float(run_value(cfg, case, "dump_interval_early_ps")), dt)
    return f"""# V8 TTM-MD periodic slab diagnostic, not strict free-surface ablation
# case: {name}
# fix ttm/mod requires fully periodic boundaries; this case verifies
# electron-lattice coupling and temperature evolution only.

{preamble}
fix             twotemp all ttm/mod {int(grid['seed'])} v80_{name}.ttm_mod {int(cfg['sizes'][case['size']]['grid_nx'])} {int(cfg['sizes'][case['size']]['grid_ny'])} {int(cfg['sizes'][case['size']]['grid_nz'])} set {float(grid['initial_electron_temperature_K']):.3f} outfile {dump_early} Te_out
{dump_and_run_lines(case, cfg, 'v80', prod_steps)}
unfix           twotemp
"""


def free_surface_laser_params(cfg: dict, case: dict) -> dict:
    laser = cfg["laser"]
    return {
        "heat_scale": float(case.get("heat_scale", laser.get("heat_scale", 1.0))),
        "n_layers": int(case.get("free_surface_heat_layers", laser["free_surface_heat_layers"])),
        "dz_A": float(case.get("free_surface_heat_layer_thickness_A", laser["free_surface_heat_layer_thickness_A"])),
        "delta_A": float(case.get("optical_penetration_depth_nm", laser["optical_penetration_depth_nm"])) * 10.0,
        "lattice_fraction": float(
            case.get("free_surface_lattice_deposition_fraction", laser["free_surface_lattice_deposition_fraction"])
        ),
        "fwhm_ps": float(laser["pulse_fwhm_ps"]),
        "center_ps": float(laser["pulse_center_ps"]),
        "layer_caps_eV": case.get("free_surface_layer_max_step_eV"),
        "default_cap_eV": float(
            case.get(
                "free_surface_max_layer_energy_per_step_eV",
                laser.get("free_surface_max_layer_energy_per_step_eV", 0.0),
            )
        ),
    }


def heat_layer_fixes(cfg: dict, case: dict, top_A: float, area_A2: float) -> str:
    lp = free_surface_laser_params(cfg, case)
    dt = case_timestep_ps(cfg, case)
    absorbed = absorbed_fluence(cfg, case)
    absorbed_eV_A2 = absorbed * J_CM2_TO_EV_A2
    total_lattice_eV = absorbed_eV_A2 * area_A2 * lp["lattice_fraction"] * lp["heat_scale"]
    n_layers = lp["n_layers"]
    dz_A = lp["dz_A"]
    delta_A = lp["delta_A"]
    sigma = gaussian_sigma_ps(lp["fwhm_ps"])
    center = lp["center_ps"]
    norm = 1.0 / (sigma * math.sqrt(2.0 * math.pi))
    layer_caps = lp["layer_caps_eV"]
    if layer_caps is None:
        default_cap = lp["default_cap_eV"]
        layer_caps = [default_cap] * n_layers if default_cap > 0 else [0.0] * n_layers
    else:
        layer_caps = list(layer_caps)
        if len(layer_caps) < n_layers:
            layer_caps.extend([layer_caps[-1] if layer_caps else 0.0] * (n_layers - len(layer_caps)))
        layer_caps = layer_caps[:n_layers]

    weights = []
    for i in range(n_layers):
        lo = i * dz_A
        hi = (i + 1) * dz_A
        weights.append(math.exp(-lo / delta_A) - math.exp(-hi / delta_A))
    wsum = sum(weights) or 1.0
    weights = [w / wsum for w in weights]

    lines = [
        "# Equivalent free-surface laser heat source.",
        "# This is not fix ttm/mod.  It deposits a calibrated fraction of absorbed",
        "# Beer-Lambert energy directly into near-surface lattice kinetic energy.",
        f"# heat_scale={lp['heat_scale']:.6g}, n_layers={n_layers}, dz_A={dz_A:.6g}, delta_A={delta_A:.6g}",
        f"# Per-layer peak-rate caps applied at build time via max_step_eV / (norm*dt).",
        f"variable        pulse equal {norm:.12g}*exp(-0.5*((time-{center:.12g})/{sigma:.12g})^2)",
    ]
    for i, weight in enumerate(weights, start=1):
        depth_lo = (i - 1) * dz_A
        depth_hi = i * dz_A
        zlo = max(0.0, top_A - depth_hi)
        zhi = max(0.0, top_A - depth_lo)
        layer_energy = total_lattice_eV * weight
        cap_i = float(layer_caps[i - 1])
        if cap_i > 0:
            cap_energy = cap_i / max(norm * dt, 1e-30)
            if layer_energy > cap_energy:
                lines.append(f"# heat_{i:02d} capped: raw={layer_energy:.6g} eV -> {cap_energy:.6g} eV (max {cap_i:.6g} eV/step)")
                layer_energy = cap_energy
        lines.extend([
            f"region          heat_{i:02d} block INF INF INF INF {zlo:.6f} {zhi:.6f} units box",
            f"group           heat_{i:02d} region heat_{i:02d}",
            f"variable        qheat_{i:02d} equal {layer_energy:.12g}*v_pulse",
            f"fix             heatfix_{i:02d} heat_{i:02d} heat 1 v_qheat_{i:02d}",
        ])
    return "\n".join(lines) + "\n"


def heat_layer_unfixes(cfg: dict, case: dict) -> str:
    n_layers = free_surface_laser_params(cfg, case)["n_layers"]
    return "\n".join(f"unfix           heatfix_{i:02d}" for i in range(1, n_layers + 1)) + "\n"


def free_surface_input(cfg: dict, case: dict) -> str:
    name = case["name"]
    dt = case_timestep_ps(cfg, case)
    prod_steps = steps(float(case["production_ps"]), dt)
    free_boundary = str(cfg["boundary"].get("free_surface_lammps_boundary", "p p fs"))
    preamble, meta = common_preamble(cfg, case, boundary=free_boundary)
    return f"""# V8 free-surface ablation morphology case
# case: {name}
# boundary {free_boundary}, no fix ttm/mod.
# z lower face is fixed; z upper face shrink-wraps so ejecta is not lost.
# Main result for surface height, crater, ejecta and plume visualization.

{preamble}
{heat_layer_fixes(cfg, case, meta['top_A'], meta['area_A2'])}
{dump_and_run_lines(case, cfg, 'v80', prod_steps)}
{heat_layer_unfixes(cfg, case)}
"""


def lammps_input(cfg: dict, case: dict) -> str:
    mode = str(case.get("case_mode", "ttm_periodic_slab_diagnostic"))
    if mode == "ttm_periodic_slab_diagnostic":
        return ttm_diagnostic_input(cfg, case)
    if mode == "free_surface_ablation":
        return free_surface_input(cfg, case)
    raise ValueError(f"Unknown V80 case_mode: {mode}")


def write_laser_diagnostics(cfg: dict, out_dir: Path, registry: list[dict[str, str]]) -> None:
    laser = cfg["laser"]
    out_dir.mkdir(parents=True, exist_ok=True)
    fwhm = float(laser["pulse_fwhm_ps"])
    center = float(laser["pulse_center_ps"])
    sigma = gaussian_sigma_ps(fwhm)
    t = np.linspace(0.0, 100.0, 1000)
    temporal = np.exp(-0.5 * ((t - center) / sigma) ** 2)
    temporal /= np.trapezoid(temporal, t)
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.plot(t, temporal)
    ax.axvline(center, color="k", ls="--", lw=0.8)
    ax.set_title("V8 laser temporal profile, 30 ps FWHM")
    ax.set_xlabel("time / ps")
    ax.set_ylabel("normalized source")
    ax.grid(True, alpha=0.25)
    fig.savefig(out_dir / "laser_temporal_profile.png", dpi=180)
    plt.close(fig)

    z = np.linspace(0.0, 80.0, 800)
    delta = float(laser["optical_penetration_depth_nm"])
    profile = np.exp(-z / delta)
    profile /= np.trapezoid(profile, z)
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.plot(z, profile)
    ax.axvline(delta, color="k", ls="--", lw=0.8, label="penetration depth")
    ax.set_title("V8 Beer-Lambert depth profile")
    ax.set_xlabel("depth below surface / nm")
    ax.set_ylabel("normalized source")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.savefig(out_dir / "laser_source_z_profile.png", dpi=180)
    plt.close(fig)

    with (out_dir / "absorbed_energy_check.csv").open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "case",
            "case_mode",
            "incident_fluence_J_cm2",
            "reflectivity_R",
            "absorbed_fluence_J_cm2",
            "absorbed_fraction",
            "pulse_fwhm_ps",
            "penetration_depth_nm",
            "I0_eV_ps_A2",
            "free_surface_lattice_deposition_fraction",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in registry:
            writer.writerow({k: row.get(k, "") for k in fieldnames})
    print(f"[OK] {out_dir / 'laser_temporal_profile.png'}")
    print(f"[OK] {out_dir / 'laser_source_z_profile.png'}")
    print(f"[OK] {out_dir / 'absorbed_energy_check.csv'}")


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--toml", type=Path, default=project_root / "config" / "v80_lammps_ttm_md_physical_model.toml")
    ap.add_argument("--out-dir", type=Path, default=project_root / "models" / "v80_lammps_ttm_md_physical_model")
    args = ap.parse_args()
    cfg = load_toml(args.toml)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    registry: list[dict[str, str]] = []
    for case in cfg["selected_cases"]:
        size = cfg["sizes"][case["size"]]
        incident = float(case["incident_fluence_J_cm2"])
        absorbed = absorbed_fluence(cfg, case)
        mode = str(case.get("case_mode", "ttm_periodic_slab_diagnostic"))
        dt = case_timestep_ps(cfg, case)
        prod_ps = float(case["production_ps"])
        run_steps = steps(prod_ps, dt)
        name = str(case["name"])
        in_path = args.out_dir / f"v80_{name}.in"
        ttm_path = args.out_dir / f"v80_{name}.ttm_mod"
        in_path.write_text(lammps_input(cfg, case), encoding="utf-8")
        if mode == "ttm_periodic_slab_diagnostic":
            ttm_path.write_text(ttm_text(cfg, case), encoding="utf-8")
            ttm_rel = str(ttm_path.relative_to(project_root)).replace("\\", "/")
            uses_ttm = "true"
            boundary = "p p p"
        else:
            if ttm_path.exists():
                ttm_path.unlink()
            ttm_rel = ""
            uses_ttm = "false"
            boundary = str(cfg["boundary"].get("free_surface_lammps_boundary", "p p fs"))

        atoms = int(size["nx_cells"]) * int(size["ny_cells"]) * int(size["nz_cells"]) * 8
        I0 = gaussian_peak_intensity(absorbed, float(cfg["laser"]["pulse_fwhm_ps"]), float(cfg["laser"]["intensity_scale"]))
        registry.append({
            "case": name,
            "name": name,
            "case_mode": mode,
            "boundary": boundary,
            "uses_fix_ttm_mod": uses_ttm,
            "size": str(case["size"]),
            "incident_fluence_J_cm2": f"{incident:.8g}",
            "reflectivity_R": f"{float(cfg['laser']['reflectivity_R']):.8g}",
            "absorbed_fluence_J_cm2": f"{absorbed:.8g}",
            "absorbed_fraction": f"{1.0 - float(cfg['laser']['reflectivity_R']):.8g}",
            "pulse_fwhm_ps": f"{float(cfg['laser']['pulse_fwhm_ps']):.8g}",
            "penetration_depth_nm": f"{float(cfg['laser']['optical_penetration_depth_nm']):.8g}",
            "I0_eV_ps_A2": f"{I0:.8g}",
            "free_surface_lattice_deposition_fraction": f"{float(cfg['laser']['free_surface_lattice_deposition_fraction']):.8g}",
            "heat_scale": f"{float(case.get('heat_scale', cfg['laser'].get('heat_scale', 1.0))):.8g}",
            "free_surface_heat_layers": str(case.get("free_surface_heat_layers", cfg["laser"]["free_surface_heat_layers"])),
            "production_ps": f"{prod_ps:.8g}",
            "timestep_ps": f"{dt:.8g}",
            "run_steps": str(run_steps),
            "expected_final_time_ps": f"{prod_ps:.8g}",
            "atoms_expected": str(atoms),
            "input_file": str(in_path.relative_to(project_root)).replace("\\", "/"),
            "ttm_mod_file": ttm_rel,
            "run_lammps_hint": str(case.get("run_lammps_hint", True)).lower(),
            "notes": str(case.get("notes", "")),
        })

    reg_path = args.out_dir / "v80_case_registry.csv"
    with reg_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(registry[0].keys()))
        writer.writeheader()
        writer.writerows(registry)
    write_laser_diagnostics(cfg, project_root / "figures" / "v80_physical_model_update" / "laser_source", registry)
    print(f"[OK] {reg_path}")
    for row in registry:
        print(
            f"[OK] {row['input_file']} mode={row['case_mode']} "
            f"boundary={row['boundary']} ttm={row['uses_fix_ttm_mod']} "
            f"steps={row['run_steps']} time={row['expected_final_time_ps']} ps"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
