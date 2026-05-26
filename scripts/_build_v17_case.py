"""
_build_v17_case.py
==================

Internal helper for V1.7 -- builds ONE LS-DYNA keyword file for the
local refined 30 ps axisymmetric thermal model.

Mesh: structured non-uniform r-z grid (NOT the V1 uniform grid).
Physics: single-temperature thermal (NOT TTM).

Consumed by scripts/build_v17_30ps_local_mesh.py.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


# =============================================================================
# Fixed-width formatters (ASCII, same as V1)
# =============================================================================
def f10(v: float) -> str:
    if v == 0.0:
        return "       0.0"
    s = f"{v:10.3e}"
    if len(s) != 10:
        s = f"{v:.2e}".rjust(10)
    return s


def i10(v: int) -> str:
    return f"{v:10d}"


# =============================================================================
# Non-uniform grid builders
# =============================================================================
def build_r_nodes_um(cfg: dict) -> np.ndarray:
    g = cfg["grid_r_um"]
    r_max = float(cfg["geometry"]["radius_um"])
    nodes = [0.0]
    r = 0.0
    while r < float(g["inner_radius_um"]) - 1e-9:
        r += float(g["inner_dr_um"])
        nodes.append(min(r, float(g["inner_radius_um"])))
    r = float(g["inner_radius_um"])
    while r < float(g["mid_radius_um"]) - 1e-9:
        r += float(g["mid_dr_um"])
        nodes.append(min(r, float(g["mid_radius_um"])))
    r = float(g["mid_radius_um"])
    while r < r_max - 1e-9:
        r += float(g["outer_dr_um"])
        nodes.append(min(r, r_max))
    nodes.append(r_max)
    return np.unique(np.round(np.asarray(nodes, dtype=float), 6))


def build_z_nodes_um(cfg: dict) -> np.ndarray:
    """
    Build z nodes from bottom (z=0) to top surface (z=thickness).
    Finest mesh in the top surface_fine_depth_um layer.
    """
    gz = cfg["grid_z_um"]
    z_top = float(cfg["geometry"]["thickness_um"])
    z_fine_start = z_top - float(gz["surface_fine_depth_um"])

    # Start from top, walk downward with prescribed dz schedule, then reverse.
    z_list = [z_top]
    z = z_top

    # Region A: top surface_fine_depth_um, constant dz
    dz_a = float(gz["surface_dz_nm"]) * 1e-3  # um
    while z > z_fine_start + 1e-9:
        z = max(z - dz_a, z_fine_start)
        z_list.append(z)

    # Region B: from z_fine_start down to (z_top - mid_depth_end_um)
    z_mid_end = z_top - float(gz["mid_depth_end_um"])
    n_mid = max(20, int((z_fine_start - z_mid_end) / 0.05))  # ~50 nm avg
    z_seg = np.linspace(z_fine_start, z_mid_end, n_mid + 1)[1:]
    z_list.extend(z_seg.tolist())

    # Region C: z_mid_end down to 0, geometrically increasing dz
    z = z_mid_end
    dz = float(gz["bulk_dz_start_nm"]) * 1e-3
    dz_max = float(gz["bulk_dz_end_nm"]) * 1e-3
    while z > 1e-9:
        dz = min(dz * 1.15, dz_max)
        z = max(z - dz, 0.0)
        z_list.append(z)

    nodes = np.unique(np.round(np.asarray(sorted(set(z_list)), dtype=float), 6))
    if nodes[0] != 0.0:
        nodes = np.concatenate([[0.0], nodes])
    if nodes[-1] != z_top:
        nodes = np.concatenate([nodes, [z_top]])
    return nodes


def mesh_stats(r_um: np.ndarray, z_um: np.ndarray) -> dict:
    dr_min = float(np.min(np.diff(r_um)))
    dr_max = float(np.max(np.diff(r_um)))
    dz_arr = np.diff(z_um)
    # minimum dz near top surface
    dz_min = float(np.min(dz_arr))
    dz_max = float(np.max(dz_arr))
    NR = len(r_um) - 1
    NZ = len(z_um) - 1
    return {
        "NR": NR, "NZ": NZ,
        "n_nodes": (NR + 1) * (NZ + 1),
        "n_elements": NR * NZ,
        "dr_min_um": dr_min, "dr_max_um": dr_max,
        "dz_min_um": dz_min, "dz_max_um": dz_max,
        "dz_min_nm": dz_min * 1000.0,
    }


# =============================================================================
# Case spec
# =============================================================================
@dataclass(frozen=True)
class V17Case:
    name: str
    Ep_uJ: float
    run_lsdyna: bool
    notes: str
    tau_ps: float
    w0_um: float
    A: float
    t_end_ns: float
    rho: float
    cp: float
    k: float
    T_init: float

    @property
    def tau_ms(self) -> float:
        return self.tau_ps * 1e-9          # ps -> ms: 30 ps = 3e-8 ms

    @property
    def w0_mm(self) -> float:
        return self.w0_um * 1e-3

    @property
    def Ep_J(self) -> float:
        return self.Ep_uJ * 1e-6

    @property
    def t_end_ms(self) -> float:
        return self.t_end_ns * 1e-6

    @property
    def I0_peak_kW_per_mm2(self) -> float:
        """
        Peak of Gaussian temporal pulse normalised to deposit Ep:

            integral_t q(t) dt = A * 2*Ep / (pi*w0^2)
            q(t) = I0 * exp(-t^2 / (2 tau^2))   (centred at t=0)
            => I0 = A * 2*Ep / (pi*w0^2 * sqrt(2*pi)*tau)
        """
        F_total = 2.0 * self.Ep_J / (math.pi * self.w0_mm ** 2)
        return self.A * F_total / (math.sqrt(2.0 * math.pi) * self.tau_ms)


@dataclass
class V17Mesh:
    r_um: np.ndarray
    z_um: np.ndarray
    stats: dict

    @property
    def r_mm(self) -> np.ndarray:
        return self.r_um * 1e-3

    @property
    def z_mm(self) -> np.ndarray:
        return self.z_um * 1e-3

    @property
    def NR(self) -> int:
        return self.stats["NR"]

    @property
    def NZ(self) -> int:
        return self.stats["NZ"]


def build_mesh(cfg: dict) -> V17Mesh:
    r_um = build_r_nodes_um(cfg)
    z_um = build_z_nodes_um(cfg)
    stats = mesh_stats(r_um, z_um)
    return V17Mesh(r_um=r_um, z_um=z_um, stats=stats)


# =============================================================================
# Node / element / flux generators
# =============================================================================
def gen_nodes(mesh: V17Mesh) -> list[str]:
    NR, NZ = mesh.NR, mesh.NZ
    out: list[str] = []
    for j in range(NZ + 1):
        y = mesh.z_mm[j]
        for i in range(NR + 1):
            x = mesh.r_mm[i]
            nid = i + j * (NR + 1) + 1
            out.append(f"{nid:8d}{x:16.7e}{y:16.7e}{0.0:16.7e}")
    return out


def gen_elements(mesh: V17Mesh, pid: int = 1) -> list[str]:
    NR, NZ = mesh.NR, mesh.NZ
    out: list[str] = []
    for j in range(NZ):
        for i in range(NR):
            eid = i + j * NR + 1
            n1 = i + j * (NR + 1) + 1
            n2 = (i + 1) + j * (NR + 1) + 1
            n3 = (i + 1) + (j + 1) * (NR + 1) + 1
            n4 = i + (j + 1) * (NR + 1) + 1
            out.append(f"{eid:8d}{pid:8d}{n1:8d}{n2:8d}{n3:8d}{n4:8d}")
    return out


def gen_top_flux_segments(mesh: V17Mesh, w0_mm: float, lcid: int) -> list[str]:
    NR, NZ = mesh.NR, mesh.NZ
    base_top = NZ * (NR + 1)
    out: list[str] = []
    CLIP = 1.0e-12
    for i in range(NR):
        r1 = mesh.r_mm[i]
        r2 = mesh.r_mm[i + 1]
        e1 = math.exp(-2.0 * (r1 / w0_mm) ** 2)
        e2 = math.exp(-2.0 * (r2 / w0_mm) ** 2)
        if e1 < CLIP and e2 < CLIP:
            continue
        m1 = -e1 if e1 >= CLIP else 0.0
        m2 = -e2 if e2 >= CLIP else 0.0
        n1 = base_top + i + 1
        n2 = base_top + i + 1 + 1
        out.append(i10(n1) + i10(n2) + i10(n2) + i10(n1))
        out.append(i10(lcid) + f10(m1) + f10(m2) + f10(m2) + f10(m1))
    return out


def gaussian_curve_points(I0: float, tau_ms: float, t_end_ms: float,
                          n_pts: int = 40) -> list[tuple[float, float]]:
    """Sample a Gaussian pulse q(t)=I0*exp(-t^2/(2 tau^2)) from t=0 to t_end."""
    pts: list[tuple[float, float]] = [(0.0, I0)]
    t_samples = np.linspace(0.0, min(5.0 * tau_ms, t_end_ms), n_pts)
    for t in t_samples[1:]:
        q = I0 * math.exp(-0.5 * (t / tau_ms) ** 2)
        pts.append((float(t), q))
    pts.append((t_end_ms, 0.0))
    # dedupe consecutive identical times
    deduped = [pts[0]]
    for t, q in pts[1:]:
        if abs(t - deduped[-1][0]) > 1e-15:
            deduped.append((t, q))
    return deduped


def resolve_database_dt_ms(cfg: dict, t_end_ms: float) -> tuple[float, str]:
    """
    Resolve *DATABASE_* output interval (ms) from config.

    LS-DYNA thermal output uses a single constant DT on:
      *DATABASE_BINARY_D3PLOT, *DATABASE_BINARY_D3THDT,
      *DATABASE_TPRINT, *DATABASE_GLSTAT

    Config [database_output].dt_output_ns (nanoseconds) is the primary control.
    Legacy fallback: solver.n_plots_desired -> t_end / n_plots.
    """
    db = cfg.get("database_output", {})
    if "dt_output_ns" in db:
        dt_ns = float(db["dt_output_ns"])
        if dt_ns <= 0:
            raise ValueError(f"database_output.dt_output_ns must be > 0, got {dt_ns}")
        strategy = str(db.get("strategy", "uniform"))
        n_frames = int(t_end_ms / (dt_ns * 1e-6)) + 1
        desc = (
            f"uniform dt={dt_ns:g} ns ({dt_ns * 1e3:g} ps), "
            f"~{n_frames} frames over {t_end_ms * 1e6:g} ns, strategy={strategy}"
        )
        return dt_ns * 1e-6, desc

    sol = cfg.get("solver", {})
    if "n_plots_desired" in sol:
        n = float(sol["n_plots_desired"])
        if n <= 0:
            raise ValueError(f"solver.n_plots_desired must be > 0, got {n}")
        dt_ms = t_end_ms / n
        dt_ns = dt_ms * 1e6
        desc = f"legacy t_end/n_plots_desired -> dt={dt_ns:.4g} ns (~{dt_ns * 1e3:.4g} ps)"
        return dt_ms, desc

    raise KeyError(
        "config must define [database_output] dt_output_ns "
        "(or legacy solver.n_plots_desired)"
    )


# =============================================================================
# .k assembly
# =============================================================================
def build_k_text(c: V17Case, mesh: V17Mesh, cfg: dict) -> str:
    nodes = gen_nodes(mesh)
    elems = gen_elements(mesh)
    flux = gen_top_flux_segments(mesh, w0_mm=c.w0_mm, lcid=1)

    I0_peak = c.I0_peak_kW_per_mm2
    sol = cfg["solver"]
    its_frac = float(sol["its_frac_of_tau"])
    its = max(c.tau_ms * its_frac, 1e-12)
    tmin = its * 1e-3
    tmax = c.tau_ms
    dt_plot, dt_desc = resolve_database_dt_ms(cfg, c.t_end_ms)

    ctrl_term = f10(c.t_end_ms)
    therm_solver = (
        i10(1) + i10(0) + i10(11) + f10(1e-6) + i10(8)
        + f10(1.0) + f10(0.0) + f10(0.0)
    )
    therm_step = (
        i10(1) + f10(0.5) + f10(its) + f10(tmin) + f10(tmax)
        + f10(50.0) + f10(0.5) + i10(0)
    )
    db_dt = f10(dt_plot)
    db_extent_card1 = "".join(i10(0) for _ in range(8))
    db_extent_card2 = "".join(i10(0) for _ in range(8))
    db_extent_card3 = (
        i10(0) + i10(0) + f10(0.0) + i10(0) + i10(0) + i10(2)
        + i10(0) + i10(0)
    )

    mat_card1 = i10(1) + f10(c.rho)
    mat_card2 = f10(c.cp) + f10(c.k)
    sec_card = i10(1) + i10(15)
    sec_card2 = f10(1.0) + f10(1.0) + f10(1.0) + f10(1.0)
    part_card = (
        i10(1) + i10(1) + i10(0) + i10(0) + i10(0)
        + i10(0) + i10(0) + i10(1)
    )
    init_temp = i10(0) + f10(c.T_init) + i10(0)

    curve_pts = gaussian_curve_points(I0_peak, c.tau_ms, c.t_end_ms)
    curve_lines = "\n".join(f"{t:20.7e}{q:20.7e}" for t, q in curve_pts)

    st = mesh.stats
    title = f"V17 {c.name} tau={c.tau_ps:g}ps Ep={c.Ep_uJ:g}uJ local mesh"[:80]

    header = f"""$ ============================================================================
$  v17_{c.name}.k
$  Generated by scripts/build_v17_30ps_local_mesh.py -- DO NOT HAND-EDIT.
$  V1.7 local refined 30 ps axisymmetric SINGLE-TEMPERATURE thermal model.
$  NOT a full TTM -- compare results with V1.6 TTM lattice Tl only.
$
$  Laser:
$      tau = {c.tau_ps:g} ps,  w0 = {c.w0_um:g} um,  A = {c.A:g},  Ep = {c.Ep_uJ:g} uJ
$      I0 (Gaussian peak) = {I0_peak:.3e} kW/mm^2
$      t_end = {c.t_end_ns:g} ns
$
$  Database output (constant DT on all *DATABASE_* cards):
$      {dt_desc}
$      DT in .k = {dt_plot:.3e} ms  (= {dt_plot * 1e6:g} ns)
$
$  Local mesh:
$      r = 0 .. {mesh.r_um[-1]:g} um,  z = 0 .. {mesh.z_um[-1]:g} um
$      NR x NZ = {st['NR']} x {st['NZ']},  nodes = {st['n_nodes']},  shells = {st['n_elements']}
$      dr = {st['dr_min_um']:.4g} .. {st['dr_max_um']:.4g} um
$      dz = {st['dz_min_nm']:.4g} .. {st['dz_max_um']*1000:.4g} nm (min at top surface)
$
$  Units: mm - ms - kg - kN - GPa - K
$ ============================================================================
*KEYWORD
*TITLE
{title}
$
*CONTROL_SOLUTION
{i10(1)}
$
*CONTROL_TERMINATION
{ctrl_term}
$
*CONTROL_THERMAL_SOLVER
{therm_solver}
$
*CONTROL_THERMAL_TIMESTEP
{therm_step}
$
$ --- database output: DT={dt_plot:.3e} ms ({dt_plot * 1e6:g} ns) ---
*DATABASE_BINARY_D3PLOT
{db_dt}
$
*DATABASE_BINARY_D3THDT
{db_dt}
$
*DATABASE_EXTENT_BINARY
{db_extent_card1}
{db_extent_card2}
{db_extent_card3}
$
*DATABASE_TPRINT
{db_dt}
$
*DATABASE_GLSTAT
{db_dt}
$
*MAT_THERMAL_ISOTROPIC
{mat_card1}
{mat_card2}
$
*SECTION_SHELL
{sec_card}
{sec_card2}
$
*PART
Silicon V17 local 30ps {c.name}
{part_card}
$
*NODE
"""

    body = (
        header
        + "\n".join(nodes) + "\n$\n"
        + "*ELEMENT_SHELL\n"
        + "\n".join(elems) + "\n$\n"
        + f"""*INITIAL_TEMPERATURE_SET
{init_temp}
$
*DEFINE_CURVE
{i10(1)}
{curve_lines}
$
*BOUNDARY_FLUX_SEGMENT
"""
        + "\n".join(flux) + "\n$\n"
        + "*END\n"
    )
    return body


def write_k(text: str, out_path: Path) -> int:
    try:
        text.encode("ascii")
    except UnicodeEncodeError as exc:
        raise RuntimeError(f"non-ASCII in .k: {text[exc.start:exc.end]!r}") from exc
    out_path.parent.mkdir(parents=True, exist_ok=True)
    body = text.replace("\r\n", "\n").replace("\n", "\r\n").encode("ascii")
    out_path.write_bytes(body)
    return len(body)


def node_id_at(mesh: V17Mesh, i: int, j: int) -> int:
    """Grid index (i=r, j=z) -> LS-DYNA node id (1-based)."""
    return i + j * (mesh.NR + 1) + 1


def axis_top_node_id(mesh: V17Mesh) -> int:
    """r=0, z=top surface node."""
    return node_id_at(mesh, 0, mesh.NZ)
