"""
内部工具：生成 models/v1_thermal_silicon.k

用法（项目 venv）：
    .\\.venv\\Scripts\\python.exe scripts\\_build_v1_mesh.py

说明：
    - 一次性运行，产出最终的 .k 文件供 LS-DYNA 直接读取
    - 单位制：mm-ms-kg-kN-GPa-K（与项目 test.k 保持一致）
    - 2D 轴对称：Y 轴为旋转轴，X 为径向，所有节点 X>=0、Z=0
    - 顶面（y=THICKNESS）施加高斯热流 + 矩形时间脉冲

可在文件顶部修改物理 / 几何参数。
"""

from __future__ import annotations

import math
from pathlib import Path

# =========================================================================
# 物理与几何参数（mm-ms-kg-kN-GPa-K 单位制）
# =========================================================================

# ---- 几何 ----
RADIUS    = 0.500       # 硅片半径   500 um = 0.5 mm
THICKNESS = 0.200       # 硅片厚度   200 um = 0.2 mm
NR = 100                # 径向单元数  -> dr = 5.0 um
NZ = 80                 # 厚度方向单元数 -> dz = 2.5 um

# ---- 硅的热物性（室温常物性，简化处理）----
RHO_SI   = 2.33e-6      # 密度        2330 kg/m^3 = 2.33e-6 kg/mm^3
CP_SI    = 700.0        # 比热容      700 J/(kg*K)
K_SI     = 1.50e-4      # 导热系数    150 W/(m*K) = 1.5e-4 kW/(mm*K)

# ---- 激光参数 ----
W0 = 0.020              # 1/e^2 光斑半径 20 um = 0.02 mm
TAU_PULSE = 1.0e-4      # 脉冲宽度    100 ns = 1e-4 ms
I0_PEAK = 10.0          # 峰值表面热流 10 kW/mm^2 = 1e10 W/m^2 （保守值）

# ---- 时间 ----
T_INIT  = 300.0         # 初始温度 K
T_END   = 5.0e-4        # 总时长 500 ns = 5e-4 ms
DT_PLOT = 2.0e-5        # d3plot 输出间隔 20 ns

# ---- 输出文件 ----
OUT = Path(__file__).resolve().parents[1] / "models" / "v1_thermal_silicon.k"


# =========================================================================
# 工具：10 列定宽数值
# =========================================================================
def f10(v: float) -> str:
    """LS-DYNA 10-char 数值字段（科学计数法，右对齐）。

    用 ``:10.3e`` 既能保证 ``+1.234e+05`` 这样的负值/正值都正好 10 字符，
    又能容纳 ``e-100`` ~ ``e+100`` 量级；零特殊处理。
    """
    if v == 0.0:
        return "       0.0"
    s = f"{v:10.3e}"  # 总宽 10：" 1.000e-06" / "-1.000e-06"
    if len(s) != 10:
        # 极端值（如 1e+100）走更紧凑写法
        s = f"{v:.2e}".rjust(10)
    return s


def i10(v: int) -> str:
    return f"{v:10d}"


# =========================================================================
# 网格生成
# =========================================================================
def gen_nodes() -> list[str]:
    dr = RADIUS / NR
    dz = THICKNESS / NZ
    out: list[str] = []
    # 节点 ID 规则：nid = i + j*(NR+1) + 1   ( i = 0..NR, j = 0..NZ )
    # x = i*dr  (径向)；y = j*dz (轴向，0=底面, NZ*dz=顶面)；z = 0
    for j in range(NZ + 1):
        y = j * dz
        for i in range(NR + 1):
            x = i * dr
            nid = i + j * (NR + 1) + 1
            out.append(f"{nid:8d}{x:16.7e}{y:16.7e}{0.0:16.7e}")
    return out


def gen_elements(pid: int = 1) -> list[str]:
    out: list[str] = []
    for j in range(NZ):
        for i in range(NR):
            eid = i + j * NR + 1
            n1 = i + j * (NR + 1) + 1                # (i,   j)
            n2 = (i + 1) + j * (NR + 1) + 1          # (i+1, j)
            n3 = (i + 1) + (j + 1) * (NR + 1) + 1    # (i+1, j+1)
            n4 = i + (j + 1) * (NR + 1) + 1          # (i,   j+1)
            out.append(f"{eid:8d}{pid:8d}{n1:8d}{n2:8d}{n3:8d}{n4:8d}")
    return out


def gen_top_flux_segments(lcid: int) -> list[str]:
    """
    顶面 (j = NZ) 上每条线段一项 *BOUNDARY_FLUX_SEGMENT。

    LS-DYNA 2D 轴对称中段定义为线段（2 节点）；为兼容 4 节点段格式，
    采用退化四边形：N1, N2, N3=N2, N4=N1，对应 mlc1..mlc4 同样退化。

    符号约定：
        总热流  = LCID(t) * mlc_i
        负值 => 热流入物体（加热）；本处用 mlc_i = -exp(-2 r^2 / w0^2)。
    """
    dr = RADIUS / NR
    base_top = NZ * (NR + 1)  # 顶面节点 ID 偏移：(0, NZ) 节点 = base_top + 1
    out: list[str] = []
    # 远离光斑中心 (r >> w0) 时 exp(-2 r^2/w0^2) 早就低于浮点下溢，
    # 但仍生成关键字会让 LS-DYNA 多解析若干段无意义的边界条件，
    # 因此把 |m| < CLIP 的值直接归零。
    CLIP = 1.0e-12
    for i in range(NR):
        r1 = i * dr
        r2 = (i + 1) * dr
        e1 = math.exp(-2.0 * (r1 / W0) ** 2)
        e2 = math.exp(-2.0 * (r2 / W0) ** 2)
        if e1 < CLIP and e2 < CLIP:
            # 整段对加热几乎无贡献，跳过
            continue
        m1 = -e1 if e1 >= CLIP else 0.0
        m2 = -e2 if e2 >= CLIP else 0.0
        n1 = base_top + i + 1
        n2 = base_top + (i + 1) + 1
        out.append(i10(n1) + i10(n2) + i10(n2) + i10(n1))
        out.append(i10(lcid) + f10(m1) + f10(m2) + f10(m2) + f10(m1))
    return out


# =========================================================================
# 拼装关键字文件
# =========================================================================
def build_k_text() -> str:
    nodes_lines = gen_nodes()
    elems_lines = gen_elements(pid=1)
    flux_lines = gen_top_flux_segments(lcid=1)

    n_nodes = len(nodes_lines)
    n_elems = len(elems_lines)
    dr_um = RADIUS / NR * 1000.0
    dz_um = THICKNESS / NZ * 1000.0

    # ---- 控制卡（10 列定宽） ----
    ctrl_term = i10(0) + f10(T_END)  # endtim only (rest defaulted)
    ctrl_term = f10(T_END)  # field 1: endtim

    therm_solver = (
        i10(1) + i10(0) + i10(11) + f10(1e-6) + i10(8)
        + f10(1.0) + f10(0.0) + f10(0.0)
    )  # atype ptype solver cgtol gpt eqheat fwork sbc

    # *CONTROL_THERMAL_TIMESTEP card 1:
    #   TS  : 0=fixed, 1=variable
    #   TIP : time-integration parameter theta (0.5 = Crank-Nicolson)
    #   ITS : initial thermal time step  (MUST be > 0)
    #   TMIN, TMAX : variable-step bounds
    #   DTEMP : max temperature change per step (variable-step driver)
    #   TSCP  : optimal CG iteration count
    #   LCTS  : load curve for time step (0 = unused)
    therm_step = (
        i10(1)       # TS  = 1, variable
        + f10(0.5)   # TIP = 0.5, Crank-Nicolson
        + f10(1e-8)  # ITS = 1e-8 ms (= 10 ns) initial step
        + f10(1e-11) # TMIN
        + f10(1e-7)  # TMAX
        + f10(50.0)  # DTEMP (max temp change per step in K)
        + f10(0.5)   # TSCP
        + i10(0)     # LCTS
    )

    db_dt = f10(DT_PLOT)

    mat_card1 = i10(1) + f10(RHO_SI)               # tmid, tro
    mat_card2 = f10(CP_SI) + f10(K_SI)             # hc, tc

    sec_card = i10(1) + i10(15)                    # secid, elform=15 (axisym vol-weighted)
    # *SECTION_SHELL 即使 ELFORM=15（轴对称体）也需要 Card 2（厚度），
    # 否则 LS-DYNA 会把下一行（*PART 的 title）当成 Card 2 读，引发连锁错误。
    # 对 ELFORM=15 来说 T1..T4 实际被忽略，给 1.0 占位即可。
    sec_card2 = f10(1.0) + f10(1.0) + f10(1.0) + f10(1.0)
    part_card = (
        i10(1) + i10(1) + i10(0) + i10(0) + i10(0)
        + i10(0) + i10(0) + i10(1)
    )  # pid secid mid eosid hgid grav adpopt tmid

    init_temp = i10(0) + f10(T_INIT) + i10(0)      # nsid=0 (all), temp, loc

    # ---- 时间脉冲（4 个折点）----
    curve_pts = [
        (0.0,                   I0_PEAK),
        (TAU_PULSE,             I0_PEAK),
        (TAU_PULSE * 1.0001,    0.0),
        (T_END,                 0.0),
    ]
    curve_lines = []
    for t, q in curve_pts:
        curve_lines.append(f"{t:20.7e}{q:20.7e}")
    curve_block = "\n".join(curve_lines)

    header = f"""$ ============================================================================
$  v1_thermal_silicon.k
$  Version 1 -- 2D axisymmetric silicon wafer, pure thermal laser heating.
$
$  Goal: validate Gaussian laser source + Si thermal properties + T field.
$        No erosion, no SPH, no plume, no structured-light post-processing.
$
$  Units: mm - ms - kg - kN - GPa - K
$     length mm   time ms   mass kg   force kN   stress GPa   temperature K
$     energy J (= kN*mm)   power kW (= kN*mm/ms)   heat flux kW/mm^2
$
$  Geometry (2D axisymmetric, Y axis = symmetry axis, X = radial):
$     radius R = {RADIUS:.3f} mm = {RADIUS*1000:.1f} um
$     thickness H = {THICKNESS:.3f} mm = {THICKNESS*1000:.1f} um
$     radial  NR = {NR} elements, dr = {dr_um:.3f} um
$     axial   NZ = {NZ} elements, dz = {dz_um:.3f} um
$     nodes  : {n_nodes}
$     shells : {n_elems}
$
$  Material (silicon, simplified room-T constant properties):
$     rho = {RHO_SI:.3e} kg/mm^3  (= 2330 kg/m^3)
$     cp  = {CP_SI:.1f} J/(kg*K)
$     k   = {K_SI:.3e} kW/(mm*K)  (= 150 W/(m*K))
$
$  Laser:
$     w0 (1/e^2 radius) = {W0*1000:.1f} um
$     pulse width tau   = {TAU_PULSE*1e6:.1f} ns (= {TAU_PULSE:.3e} ms)
$     peak surface flux I0 = {I0_PEAK:.2f} kW/mm^2 (= 1e10 W/m^2, conservative)
$     spatial : q(r) = I0 * exp(-2 r^2 / w0^2)
$     temporal: rectangular pulse on [0, tau]
$
$  Sign convention for *BOUNDARY_FLUX in LS-DYNA: negative value => heat
$  flows INTO the body (heating). Here mlc_i = -exp(-2 r^2 / w0^2) combined
$  with a positive LCID curve produces heating.
$ ============================================================================
*KEYWORD
*TITLE
V1 2D Axisymmetric Silicon Laser Thermal
$
$ ---------------------- Solution control ----------------------
*CONTROL_SOLUTION
$#    soln       nlq     isnan     lcint     lcacc
{i10(1)}
$
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol
{ctrl_term}
$
$ ---------------------- Thermal solver ----------------------
*CONTROL_THERMAL_SOLVER
$#   atype     ptype    solver     cgtol       gpt    eqheat     fwork       sbc
{therm_solver}
$
*CONTROL_THERMAL_TIMESTEP
$#      ts       tip       its      tmin      tmax     dtemp      tscp      lcts
{therm_step}
$
$ ---------------------- Output (binary + ASCII) ----------------------
*DATABASE_BINARY_D3PLOT
$#      dt
{db_dt}
$
*DATABASE_BINARY_D3THDT
$#      dt
{db_dt}
$
*DATABASE_TPRINT
$#      dt
{db_dt}
$
*DATABASE_GLSTAT
$#      dt
{db_dt}
$
$ ---------------------- Thermal material (silicon, constant) ----------------------
*MAT_THERMAL_ISOTROPIC
$#    tmid       tro     tgrlc    tgmult      tlat      hlat
{mat_card1}
$#      hc        tc
{mat_card2}
$
$ ---------------------- Section / Part ----------------------
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt    qr/irid     icomp     setyp
{sec_card}
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
{sec_card2}
$
*PART
Silicon wafer 2D axisymmetric
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
{part_card}
$
$ ---------------------- Nodes ----------------------
*NODE
"""

    body_nodes = "\n".join(nodes_lines) + "\n$\n"

    elem_header = "*ELEMENT_SHELL\n"
    body_elems = "\n".join(elems_lines) + "\n$\n"

    init_block = f"""$ ---------------------- Initial temperature (all nodes) ----------------------
*INITIAL_TEMPERATURE_SET
$#    nsid      temp       loc
{init_temp}
$
"""

    flux_block = f"""$ ---------------------- Laser temporal pulse: rectangular [0, tau] ----------------------
$  LCID 1 = peak surface heat flux vs time (kW/mm^2).
$  Multiplied by per-segment mlc = -exp(-2 r^2 / w0^2) gives the Gaussian heating.
*DEFINE_CURVE
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp
{i10(1)}
$#               a1                  o1
{curve_block}
$
$ ---------------------- Top-surface Gaussian heat flux (per-segment) ----------------------
*BOUNDARY_FLUX_SEGMENT
$#      n1        n2        n3        n4
$#    lcid      mlc1      mlc2      mlc3      mlc4
"""
    body_flux = "\n".join(flux_lines) + "\n$\n"

    footer = "*END\n"

    return (
        header
        + body_nodes
        + elem_header
        + body_elems
        + init_block
        + flux_block
        + body_flux
        + footer
    )


def main() -> None:
    txt = build_k_text()
    # LS-DYNA 的解析器对非 ASCII 字符（甚至 $ 注释里的中文）极其不友好，
    # 实测会触发 "input error found in structured input"，所以强制写 ASCII。
    try:
        txt.encode("ascii")
    except UnicodeEncodeError as exc:
        offending_line = txt[: exc.start].count("\n") + 1
        raise RuntimeError(
            f"非 ASCII 字符在第 {offending_line} 行: {txt[exc.start:exc.end]!r}.\n"
            "请把它改成纯 ASCII 后再生成 .k 文件。"
        ) from exc

    OUT.parent.mkdir(parents=True, exist_ok=True)
    # 写为带 CRLF 的 ASCII，避免 Windows 端混合行尾干扰
    OUT.write_bytes(txt.replace("\r\n", "\n").replace("\n", "\r\n").encode("ascii"))
    n_nodes = (NR + 1) * (NZ + 1)
    n_elems = NR * NZ
    size_kb = OUT.stat().st_size / 1024
    print(f"Wrote {OUT}")
    print(f"  nodes    : {n_nodes}")
    print(f"  elements : {n_elems}")
    print(f"  size     : {size_kb:.1f} KB")


if __name__ == "__main__":
    main()
