"""
check_v1_outputs.py
====================

检查 Version 1（2D 轴对称硅片纯热模型）的运行结果是否合理。

用法（项目 venv，项目根目录下运行）：
    .\\.venv\\Scripts\\python.exe scripts\\check_v1_outputs.py

默认结果目录：results\\v1\\
也可以传入自定义目录：
    .\\.venv\\Scripts\\python.exe scripts\\check_v1_outputs.py results\\v1

本脚本不依赖二进制 d3plot 解析（避免引入额外库）。它做的事：

[A] 文件存在性
    - d3plot          → 二进制温度场
    - d3hsp           → 求解器回显（确认参数被正确吃进）
    - messag          → 警告 / 错误
    - glstat          → 全局统计（每个时间步）
    - tprint          → 节点 / 单元温度 ASCII 时程
    - bndout          → 边界条件输出（含热流）— 可选

[B] 正常终止
    - messag 末尾出现 "Normal termination"

[C] 求解器健康
    - messag / d3hsp 中 ERROR / FATAL 计数
    - WARNING 计数（可容忍但要看）

[D] 关键参数回显
    - 在 d3hsp 中找回我们写进 .k 的几个数字
      （材料密度、比热、热导率、激光峰值热流、脉宽…）

[E] 解析估计 vs 物理直觉
    - 给出"理论峰值温升 ΔT" 与"室温起点 300 K 后的预期峰值温度"
    - 用户从 LS-PrePost / d3plot 读取的真实峰值温度可以对照本估计

本脚本 NOT 做的事（V1 不需要）：
    - 不读取 d3plot 二进制
    - 不画温度云图
    - 不算梯度
"""

from __future__ import annotations

import math
import os
import re
import sys
from pathlib import Path

# --------- Windows 控制台 UTF-8 -----------------------------------------------
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass
    try:
        os.system("chcp 65001 > nul")
    except Exception:
        pass


# =============================================================================
# 必须与 v1_thermal_silicon.k 保持一致（仅用于解析估计 + 回显比对）
# =============================================================================
PARAMS = {
    "unit_system": "mm-ms-kg-kN-GPa-K",
    # 几何
    "R_mm": 0.500,
    "H_mm": 0.200,
    "NR": 100,
    "NZ": 80,
    # 材料 — 硅 (室温常物性)
    "rho_kg_per_mm3": 2.33e-6,
    "rho_SI_kg_per_m3": 2330.0,
    "cp_J_per_kgK": 700.0,
    "k_kW_per_mmK": 1.50e-4,
    "k_SI_W_per_mK": 150.0,
    # 激光
    "w0_mm": 0.020,
    "tau_pulse_ms": 1.0e-4,
    "tau_pulse_s": 1.0e-7,
    "I0_kW_per_mm2": 10.0,
    "I0_SI_W_per_m2": 1.0e10,
    # 时间
    "T_init_K": 300.0,
    "T_end_ms": 5.0e-4,
}

# 硅的熔点 / 沸点（参考）
T_MELT_SI_K = 1687.0
T_BOIL_SI_K = 3538.0

SEP = "=" * 72
RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    flag = "[ OK ]" if ok else "[FAIL]"
    print(f"  {flag}  {name}")
    if detail:
        for line in detail.splitlines():
            print(f"          {line}")


def section(title: str) -> None:
    print(f"\n{SEP}\n  {title}\n{SEP}")


# =============================================================================
# [A] 文件存在性
# =============================================================================
EXPECTED_FILES = {
    "d3plot":  ("二进制温度场",   True),
    "d3hsp":   ("求解器回显",     True),
    "messag":  ("警告 / 错误",    True),
    "glstat":  ("全局统计",       False),  # GLSTAT 在纯热问题里有时不写
    "tprint":  ("温度 ASCII 输出", False),
    "bndout":  ("边界条件输出",   False),
}


def check_files(out_dir: Path) -> None:
    section(f"[A] 检查输出文件 (在 {out_dir})")
    for name, (desc, required) in EXPECTED_FILES.items():
        p = out_dir / name
        if p.is_file():
            kb = p.stat().st_size / 1024
            record(f"{name:8s}  ({desc})", True, f"size = {kb:.1f} KB")
        else:
            note = "缺失"
            if not required:
                note += "（非强制）"
            record(f"{name:8s}  ({desc})", required is False, note)


# =============================================================================
# [B] Normal termination
# =============================================================================
def check_normal_termination(messag: Path) -> None:
    section("[B] Normal termination")
    if not messag.is_file():
        record("messag 文件存在", False, f"未找到 {messag}")
        return
    txt = messag.read_text(encoding="utf-8", errors="ignore")
    # LS-DYNA 用带空格的"展开"格式打印横幅：
    #   ' N o r m a l    t e r m i n a t i o n'
    # 先把所有空白压掉再做大小写不敏感的匹配。
    compact = re.sub(r"\s+", "", txt).lower()
    ok = "normaltermination" in compact
    record("messag 中找到 'Normal termination'", ok,
           "未找到 → 检查 messag 末尾，可能仍在跑或异常退出" if not ok else "")


# =============================================================================
# [C] 求解器健康（ERROR / WARNING 计数）
# =============================================================================
def check_solver_health(out_dir: Path) -> None:
    section("[C] 求解器健康检查")
    for fname in ("messag", "d3hsp"):
        p = out_dir / fname
        if not p.is_file():
            record(f"{fname} 健康检查", False, "文件不存在")
            continue
        txt = p.read_text(encoding="utf-8", errors="ignore")
        n_err = len(re.findall(r"\bERROR\b", txt))
        n_fat = len(re.findall(r"\bFATAL\b", txt, flags=re.IGNORECASE))
        n_warn = len(re.findall(r"\bWARNING\b", txt))
        ok = (n_err == 0 and n_fat == 0)
        record(
            f"{fname}: ERROR={n_err}  FATAL={n_fat}  WARNING={n_warn}",
            ok,
            "" if ok else "存在 ERROR/FATAL，需要看具体内容！",
        )


# =============================================================================
# [D] 关键参数回显：在 d3hsp 中找几个我们写进去的数字
# =============================================================================
PARAM_PATTERNS = [
    # (描述, 正则) — 这些都是 d3hsp 真实出现过的字符串
    ("solution type = thermal only",
        re.compile(r"thermal\s+only", re.I)),
    ("element formulation 15 (axisymmetric)",
        re.compile(r"axisym|elform\s*[:=]?\s*15", re.I)),
    ("Thermal property type 1 (isotropic)",
        re.compile(r"thermal\s+property\s+type", re.I)),
    ("material density = 2.33e-6",
        re.compile(r"density\s*=\s*2\.3\d{1,3}E-0?6", re.I)),
    ("heat capacity = 700",
        re.compile(r"heat\s+capacity\s*=\s*7\.0{1,5}E\+0?2", re.I)),
    ("thermal conductivity = 1.5e-4",
        re.compile(r"thermal\s+conductivity\s*=\s*1\.50*E-0?4", re.I)),
]


def check_d3hsp_echo(out_dir: Path) -> None:
    section("[D] d3hsp 关键参数回显")
    p = out_dir / "d3hsp"
    if not p.is_file():
        record("d3hsp 存在", False, "无法做参数回显检查")
        return
    txt = p.read_text(encoding="utf-8", errors="ignore")
    for desc, pat in PARAM_PATTERNS:
        ok = bool(pat.search(txt))
        record(desc, ok, "" if ok else "在 d3hsp 中未找到（可能换了写法，不一定是错）")


# =============================================================================
# [E0] 从 tprint 提取 T_max(t) 时程
# =============================================================================
def parse_tprint_tmax(tprint: Path) -> list[tuple[float, float, int]]:
    """
    解析 ASCII tprint，返回 [(t_ms, T_max_K, node_id_of_max), ...]。

    LS-DYNA 在 tprint 中对每个输出时刻都会先打印一段：
        time =  1.0005E-04     time step = ...
        minimum temperature =  300.00000   at node     ...
        maximum temperature =  501.51419   at node     8081
    我们只取 maximum 行。
    """
    if not tprint.is_file():
        return []
    txt = tprint.read_text(encoding="utf-8", errors="ignore")
    times = re.findall(r"time\s*=\s*(\S+)\s+time\s+step", txt)
    maxes = re.findall(r"maximum\s+temperature\s*=\s*(\S+)\s+at\s+node\s+(\d+)", txt)
    n = min(len(times), len(maxes))
    out: list[tuple[float, float, int]] = []
    for i in range(n):
        try:
            t = float(times[i])
            T = float(maxes[i][0])
            nid = int(maxes[i][1])
            out.append((t, T, nid))
        except ValueError:
            pass
    return out


def print_temperature_history(tprint: Path) -> None:
    section("[E0] 温度时程（tprint 解析）")
    hist = parse_tprint_tmax(tprint)
    if not hist:
        record("tprint 解析", False, f"未能解析 {tprint}（可能没生成或格式不同）")
        return

    print(f"  {'t (ms)':>12s}   {'T_max (K)':>10s}   {'ΔT (K)':>8s}   node")
    print(f"  {'-'*12}   {'-'*10}   {'-'*8}   {'-'*5}")
    T0 = PARAMS["T_init_K"]
    peak_T = T0
    peak_t = 0.0
    peak_node = -1
    for t, T, nid in hist:
        print(f"  {t:12.4e}   {T:10.2f}   {T-T0:8.2f}   {nid}")
        if T > peak_T:
            peak_T = T
            peak_t = t
            peak_node = nid

    print()
    print(f"  全局峰值：T = {peak_T:.2f} K （ΔT = {peak_T - T0:.2f} K）"
          f"  @ t = {peak_t*1e6:.1f} ns,  node {peak_node}")
    record("从 tprint 成功提取 T_max(t)", True,
           f"采样 {len(hist)} 个时刻，全局峰值 {peak_T:.1f} K")


# =============================================================================
# [E] 解析估计：高斯激光 + 半无限体表面温升
# =============================================================================
def analytical_estimate() -> dict[str, float]:
    """
    半无限体 + 表面常值矩形热流（局部估计中心 r=0 的峰值温升）:
        ΔT_max ≈ 2 * I0 * sqrt(α * τ / π) / k
    其中：
        α  = k / (ρ * cp)   热扩散率
        I0 = 表面峰值热流（W/m^2）
        τ  = 脉冲宽度（s）
        k  = 导热系数（W/(m K)）
    """
    k   = PARAMS["k_SI_W_per_mK"]            # W/(m K)
    rho = PARAMS["rho_SI_kg_per_m3"]         # kg/m^3
    cp  = PARAMS["cp_J_per_kgK"]             # J/(kg K)
    I0  = PARAMS["I0_SI_W_per_m2"]           # W/m^2
    tau = PARAMS["tau_pulse_s"]              # s
    w0  = PARAMS["w0_mm"] * 1e-3             # m
    T0  = PARAMS["T_init_K"]

    alpha = k / (rho * cp)
    L_th = math.sqrt(alpha * tau)             # 热扩散长度 (m)
    dT_max = 2.0 * I0 * math.sqrt(alpha * tau / math.pi) / k

    # 高斯吸收的总能量
    E_abs = I0 * tau * math.pi * w0 * w0 / 2.0  # J

    return {
        "alpha_m2_per_s": alpha,
        "L_th_um": L_th * 1e6,
        "dT_max_K": dT_max,
        "T_peak_pred_K": T0 + dT_max,
        "E_abs_uJ": E_abs * 1e6,
    }


def print_analytical_estimate() -> None:
    section("[E] 解析估计（半无限体 + 矩形脉冲，仅用于直觉对照）")
    est = analytical_estimate()
    print(f"  α (热扩散率)              = {est['alpha_m2_per_s']:.3e} m^2/s")
    print(f"  脉冲期间热扩散长度 L_th    = {est['L_th_um']:.2f} um")
    print(f"  吸收总能量 E_abs           = {est['E_abs_uJ']:.3f} uJ")
    print(f"  预期峰值温升 ΔT_max       ≈ {est['dT_max_K']:.1f} K")
    print(f"  预期峰值温度 T_peak       ≈ {est['T_peak_pred_K']:.1f} K"
          f"   (初温 {PARAMS['T_init_K']:.0f} K)")

    print()
    print(f"  参考：硅熔点 ≈ {T_MELT_SI_K:.0f} K，沸点 ≈ {T_BOIL_SI_K:.0f} K")

    # 合理性提示
    Tp = est["T_peak_pred_K"]
    if Tp < 400:
        verdict = "温升偏小（< 100 K），如要更明显热效应可适当增大 I0_PEAK。"
        ok = True
    elif Tp < T_MELT_SI_K - 200:
        verdict = "理论峰值远低于熔点，V1 数值上应当稳定。"
        ok = True
    elif Tp < T_MELT_SI_K:
        verdict = "理论峰值接近熔点，仍未达到。OK 但已经偏激进。"
        ok = True
    elif Tp < T_BOIL_SI_K:
        verdict = "理论峰值已超过熔点 — V1 不带相变/侵蚀，结果会出现非物理高温。建议先调小 I0_PEAK。"
        ok = False
    else:
        verdict = "理论峰值超过沸点 — 这种情况 V1 必定数值离谱。请先减小 I0_PEAK。"
        ok = False

    record("理论峰值在合理范围（远低于硅熔点 1687 K）", ok, verdict)


# =============================================================================
# [F] 调参建议
# =============================================================================
def print_tuning_advice() -> None:
    section("[F] 如果温度偏离预期")
    print("""\
  如果 LS-PrePost 看到的峰值温度 ≪ 解析估计：
      - 检查 *BOUNDARY_FLUX_SEGMENT 的 mlc 符号是不是反了
        （V1 用 负 mlc + 正 LCID = 加热；反过来变成冷却）
      - 检查 *INITIAL_TEMPERATURE_SET 是否被识别 (d3hsp 里能找到 300 K)
      - 检查 *MAT_THERMAL_ISOTROPIC 单位是否匹配（最常见的"差几个数量级"问题）

  如果峰值温度 ≫ 解析估计（或 5xx K 远高于 800 K 预期）：
      - 减小 I0_PEAK（先调 5→2 kW/mm^2）
      - 检查脉宽 TAU_PULSE 是否被 LCID 正确表达
      - 检查网格是否过粗（顶面 dz 应远小于 L_th）

  如果想故意提高温度（V2 准备阶段）：
      - 把 I0_PEAK 改到 30–50 kW/mm^2 看是否能逼近熔点
      - 千万先把 *DAMPING* / *MAT_ADD_EROSION* 准备好（V1 不要做）
""")


# =============================================================================
# 主流程
# =============================================================================
def main(argv: list[str]) -> int:
    # 默认结果目录：项目根/results/v1
    # 用 __file__ 推断项目根，而不是 Path("results/v1") 这种依赖 cwd 的相对路径，
    # 否则在 scripts/ 下执行就会去找 scripts/results/v1，导致 [A] 全 FAIL。
    if len(argv) >= 2:
        out_dir = Path(argv[1]).resolve()
    else:
        project_root = Path(__file__).resolve().parents[1]
        out_dir = (project_root / "results" / "v1").resolve()

    print(SEP)
    print("  Version 1 — 2D 轴对称硅片纯热模型 · 输出检查")
    print(f"  结果目录 : {out_dir}")
    print(SEP)

    check_files(out_dir)
    check_normal_termination(out_dir / "messag")
    check_solver_health(out_dir)
    check_d3hsp_echo(out_dir)
    print_temperature_history(out_dir / "tprint")
    print_analytical_estimate()
    print_tuning_advice()

    section("汇总")
    n_ok = sum(1 for _, ok, _ in RESULTS if ok)
    n_total = len(RESULTS)
    print(f"  通过 {n_ok} / {n_total}")
    failed = [name for name, ok, _ in RESULTS if not ok]
    if failed:
        print("  未通过项：")
        for name in failed:
            print(f"    - {name}")
        return 1
    print("  全部检查通过 — 请进一步用 LS-PrePost 打开 d3plot 验证温度场形状。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
