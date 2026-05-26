# LS-DYNA + Python 结构光辅助硅片激光烧蚀仿真总方案

> 版本：2026-05-26 起草，与 V1 验证完成同步。
> 文档用途：作为后续 V1.5 → V2 → V3 → V4 各阶段的总体规划，
> 给后续 Cursor Agent / 协作者一个统一的物理与工程目标。

---

## 1. 项目目标

本项目的**最终交付**不是单一的 LS-DYNA 模型，而是一整套
"实验 + 仿真 + 后处理 + 可视化" 的工具链，用于支撑：

1. **硅片激光烧蚀**  
   通过强紫外 / 可见光泵浦在 intrinsic silicon 表面产生**熔融区 + 蒸发区 + 坑形**。
2. **Pump-probe 成像**  
   用第二束光（探测光，建议 800 nm 或可见光）在不同延时下捕获烧蚀过程，对应实验上的高速相机 / 条纹相机帧。
3. **Streak camera / CUP 观测**  
   将 streak / CUP 测得的强度时序与仿真预测的 plume / shock 时序对照（**Compressed Ultrafast Photography**）。
4. **结构光辅助诊断**  
   通过在探测光路上插入条纹（如正弦光栅 / DMD）做**结构光照明**，
   把热致表面位移、坑形深度、plume 折射率扰动等信息映射到**条纹相位失真 + 调制度衰减**上。
5. **最终输出四层图**（见 §8）：
   - 第 1 层：普通照明序列（强度）
   - 第 2 层：结构光原始条纹序列
   - 第 3 层：解调相位 / 调制度 / 质量图
   - 第 4 层：定量曲线（坑深、坑半径、plume 前沿、velocity、SNR/CNR）

**LS-DYNA 在这条链里的角色**：
提供 *物理一致* 的温度场 / 坑形场 / 应力波场 / 喷射场，作为驱动后续光学正向仿真的"地面真值"。
LS-DYNA **不**生成光学图像；光学部分全部由 Python 完成。

---

## 2. 当前 V1 结果总结

V1 (`models/v1_thermal_silicon.k` + `results/v1/`) 是**纯热模型**，已 Normal termination，
17/17 自动检查通过，详见 `README_v1.md`。

| 项 | 值 |
| --- | --- |
| 模型类型 | 2D 轴对称纯热传导（SOLN = 1） |
| 几何 | 硅片半径 500 µm × 厚度 200 µm |
| 网格 | 径向 100 × 轴向 80 = 8000 单元，8181 节点 |
| 激光 | 高斯 w₀ = 20 µm，矩形脉冲 τ = 100 ns，I₀ = 10 kW/mm² |
| 实测峰值温度 | **501.5 K @ t = 100.1 ns**（中心顶面节点） |
| 解析估计 | 528 K（半无限体 + 矩形脉冲），偏差 ~5% |

关键事实，给后续读者：

- V1 是 **2D 轴对称**模型，`r = 0` 位于左边界，所以**最高温出现在左上角**是
  *物理正确* 的，不是偏心。`scripts/plot_v1_full_section.py` 通过镜像
  `r > 0` 到 `r < 0` 来把它显示成"完整截面 + 中心热斑"，论文里更直观。
  `scripts/plot_v1_revolved_3d.py` 则进一步把 `T(r, z)` 绕 `z` 轴旋转，
  生成 3D 圆盘视图，把"左上角热斑"翻译成"圆盘中心热点"。
- V1 验证的是 **热源 + 单位制 + 热扩散 + 输出通道** 这一条链，
  *不*代表真实烧蚀场景。
- V1 的峰值温度 ~501 K **远低于** 硅熔点 1687 K 和沸点 3538 K，
  所以 V1 内部不会出现熔融 / 汽化，也就不需要相变 / 侵蚀 / SPH。
- V1 的 `d3plot` 已经能在 LS-PrePost 中显示 Temperature（`THERM = 2`）。

> 后续所有版本以 V1 为基线，**不再修改 V1 物理参数**（密度 / 比热 / 导热 / 几何 / 单位制）。

---

## 3. 实验参数整理（基线）

实验里关心的是**单脉冲烧蚀**而非平均功率，所以参数表用 `(Ep, fluence, A, w0, τ, λ)`
描述，把重复频率单独列出来。

### 3.1 样品

| 参数 | 值 / 范围 | 备注 |
| --- | --- | --- |
| 材料 | intrinsic silicon | 室温电阻率 ≥ 200 Ω·cm |
| 取向 | (100) 或 (111) | V1 当作各向同性 |
| 厚度 | 200 µm | 与 V1 保持一致 |
| 表面 | 化学抛光 + RCA 清洗 | 用于计算反射率 R |
| 室温物性 | ρ = 2330 kg/m³, cp = 700 J/(kg·K), k = 150 W/(m·K) | 后续 V1.5 起加 ρ(T)、cp(T)、k(T) |
| 熔点 T_melt | 1687 K | V2 阈值 |
| 沸点 T_vap | 3538 K | V2 阈值 |
| 蒸发焓 ΔH_v | 1.79 × 10⁷ J/kg (= 17.9 µJ/ng) | V3 用于估 ablation depth |

### 3.2 泵浦光（候选）

| 波长 λ | 用途 | 反射率 R (Si, 法向)* | 吸收深度 α⁻¹ |
| --- | --- | --- | --- |
| 266 nm | 强吸收，最浅 | ~0.61 | ~ 5 nm |
| 355 nm | 平衡选项 | ~0.58 | ~10 nm |
| 532 nm | 商用 DPSS | ~0.37 | ~1 µm |
| 1064 nm | 备选 | ~0.32 | ~100 µm（接近透明） |

*室温 + 抛光面、强吸收波段都是面式吸收，可继续用 V1 的 surface flux 模型。
1064 nm 接近 Si 带边，吸收深度变大、必须切到体吸收 → 后续单独建模。

吸收率 A = 1 − R，例如 266 nm → A ≈ 0.39；532 nm → A ≈ 0.63；1064 nm → A ≈ 0.68。

### 3.3 探测光（pump-probe / structured-light）

| 参数 | 推荐 | 备注 |
| --- | --- | --- |
| 波长 | 800 nm (fs/ps) 或 532 nm / 632.8 nm (CW + gated camera) | 不影响 LS-DYNA |
| 模式 | 准直平面波 + DMD 条纹掩模 | 仿真侧建模为正弦灰度 |
| 条纹频率 | 50–200 lp/mm（视成像放大率） | 用 §7 中 `f` 表征 |
| 偏振 | 线偏振 + 起偏 / 检偏 | 仿真里只走标量场 |

### 3.4 脉冲参数（扫描区间）

| 参数 | V1 基线 | V1.5 扫描区间 | 物理意义 |
| --- | --- | --- | --- |
| Ep（单脉冲能量） | — | 0.1 µJ ↔ 100 µJ | 决定总沉积能量 |
| fluence F = Ep / (π w₀²) | 7.96 J/cm²（推算） | 0.1 ↔ 100 J/cm² | 真正的"强度"指标 |
| w₀（1/e² 半径） | 20 µm | 10 ↔ 50 µm | 决定径向梯度 |
| τ（脉宽） | 100 ns | 100 ns → 30 ps | 决定热扩散是否"赶得上" |
| A（吸收率） | 1.0（V1 简化） | 0.3 ↔ 0.8（按 λ） | 表层入射→吸收的转换效率 |
| 重复频率 frep | — | 不在单脉冲扫描内 | 影响平均热积累，单发不变 |

实验里**重复频率 10 Hz 只影响平均热积累**（脉冲间冷却 ≫ 单脉冲时间），所以
单发烧蚀仿真里 frep 不进 LS-DYNA，只在多脉冲版本里加上"两脉冲间冷却 100 ms"。

### 3.5 实验现场推荐光斑直径

- 直径 70 µm → 半径 35 µm。
- V1 用 w₀ = 20 µm，是更保守的选择；V1.5 起扫描到 35 µm 以匹配实际实验。

---

## 4. V1.5：真实实验参数扫描（**下一步**）

**目的**：在不引入侵蚀 / 相变的前提下，把 V1 升级成**参数化的热模型**，
扫描 `(Ep, w0, τ, A, λ)`，找出使硅表面达到熔点 / 汽化温度的边界。

### 4.1 改动点（相对 V1）

1. 把 `_build_v1_mesh.py` 的输入参数抽到 `config/laser_cases.yaml`；
2. 每个 case 生成一份独立的 `.k`，跑独立结果目录 `results/v1p5/case_XX/`；
3. 仍用同一套 mesh / 材料；
4. **不**加侵蚀，**不**加 SPH，**不**加相变。

### 4.2 扫描矩阵（建议）

| 维度 | 取值 | # |
| --- | --- | --- |
| τ (ns) | 100, 30, 10, 3, 1, 0.3, 0.1 | 7 |
| w₀ (µm) | 10, 20, 35, 50 | 4 |
| Ep (µJ) | 0.5, 5, 50, 200 | 4 |
| A | 0.4 (266 nm), 0.6 (532 nm) | 2 |

共 7 × 4 × 4 × 2 = **224 case**，分批跑。

> 注：脉宽 < 1 ns 时热扩散长度 < 100 nm，必须**加密表层网格**
> （e.g. 顶面 2 µm 内 dz = 20 nm，过渡到内部 dz = 2.5 µm）；
> 此时网格成本暴涨，应单独做 "V1.5b 小尺度"。

### 4.3 输出 / 后处理

每个 case：

- `T_max_history.csv`：t, T_max(t), r_at_max, z_at_max
- `T_top_radial.csv`：顶面径向温度曲线（每个 case 的峰值时刻）
- `melt_threshold_hit`：bool，是否触及 1687 K
- `vapor_threshold_hit`：bool，是否触及 3538 K
- `t_first_melt`, `t_first_vapor`：首次越界时刻

汇总图：

1. `Tmax_vs_Ep.png`：固定 (τ, w₀, A)，扫 Ep；
2. `melt_threshold_map.png`：在 (Ep, τ) 平面画 "是否熔" 的边界；
3. `vapor_threshold_map.png`：同上，沸点边界；
4. `crater_proxy_radius.png`：达到熔点的最外侧半径 → 等效坑半径估计；
5. `crater_proxy_depth.png`：达到熔点的最深 z → 等效坑深估计（V2 起替换为真坑深）。

### 4.4 已知风险

- τ → 0.1 ns 时 Crank-Nicolson 时间步可能不够；要么用更小 ITS，要么切到隐式自适应；
- 大 fluence 会让中心温度看似 > 10000 K，但 V1.5 不做相变，结果**只可作阈值估计**，
  不能用于真实坑深。

---

## 5. V2：温度阈值烧蚀坑模型

**目的**：从 V1.5 的温度场推导**等效坑形** crater(t)，作为 V3 / V4 的几何输入。

### 5.1 物理模型（先简化）

不在 LS-DYNA 里做 `*MAT_ADD_EROSION` 或 SPH，
而是在 Python 后处理里用**温度阈值法**：

```
crater_mask(r, z, t) = T(r, z, t) >= T_vap        # vapor / ablation
melt_mask  (r, z, t) = T_melt <= T(r, z, t) < T_vap
```

进一步**累积**：
- 一旦某个单元被标记为 `crater`，认为它已经被剥离，下一帧仍然记为"已剥离"；
- 熔融区是可逆的，跟随温度回落。

### 5.2 输出

| 文件 | 含义 |
| --- | --- |
| `crater_depth_t.csv` | 中心轴上的坑深 d(t) |
| `crater_radius_t.csv` | 顶面上的坑径 R(t) |
| `melt_depth_t.csv` | 熔融区最深 z |
| `ablation_depth_t.csv` | 与 crater_depth 同义，单列 |
| `removed_region_map_{t}.png` | 在 r-z 截面上画 crater_mask ∪ melt_mask |
| `final_crater_section.png` | t = t_end 时的 mask（"最终坑形"） |

### 5.3 选项 B（备用，若用户后期同意做侵蚀）

`*MAT_ADD_EROSION` + `*MAT_THERMAL_ISOTROPIC_TD`（温度依赖）+
`failure on T > 3538 K`。需要再讨论稳定性，先**不**做。

### 5.4 验证标准

V2 通过的标志：
- crater_depth 单调非降；
- crater_radius 在脉冲结束后达到平台并不再扩张；
- ablation volume 与 `Ep × A / ΔH_v / ρ` 在 10× 量级内对得上（粗略守恒）。

---

## 6. V3：plume / shock / ejecta 简化动力学

V2 给出坑形，V3 给出**喷出物 + 冲击波** 的等效几何与时序。

### 6.1 路线 A（**先做**）：LS-DYNA 只做坑形，Python 后处理生成 plume

把 plume 当成一个**自相似膨胀模型**：

```
plume_front_R(t) = R0 + v0 * (t - t0) + 0.5 * a * (t - t0)^2
                                  (Sedov-like, with empirical a)
plume_density_proxy(r, z, t)
    = M_total / V(t) * shape_factor(r, z, R_front(t))
```

输入：V2 的 ablation_volume(t)、ablation_velocity(t)（由 mass-flux 反推）。
输出：`plume_front_R.csv`, `plume_velocity.csv`, `plume_density_proxy_map.npz`。

**优点**：完全 Python，跟 LS-DYNA 完全解耦，便于快速迭代。

### 6.2 路线 B（**后做**，需要更多人时）：表层 SPH + recoil pressure

仍以 V2 的坑形为底，**只在 crater_mask** 的区域生成 SPH 粒子，
赋初速 = recoil pressure 反推；在 LS-DYNA 里跑 100 ns ~ 5 µs；
输出粒子位置、速度、温度时序。

约束：

- 仅做单脉冲；
- SPH 粒子尺寸 ~ 2 µm；
- 边界用 `*BOUNDARY_NON_REFLECTING_2D` 防反射。

### 6.3 输出（路线 A）

| 文件 | 内容 |
| --- | --- |
| `plume_front_R_t.csv` | 时间 vs plume 前沿半径 |
| `plume_velocity_t.csv` | 前沿速度 |
| `shock_front_R_t.csv` | 简化激波前沿（Sedov-Taylor 估计） |
| `ejecta_mask_{t}.npz` | r-z 截面上的 ejecta mask（boolean） |
| `intensity_decay_proxy.csv` | ∫ density * absorption dr（probe 强度衰减估计） |

> 即使路线 A 不用 LS-DYNA 跑新模型，**仍然在 docs/README 里写明** "如何切到路线 B"，
> 保证后续可以无痛升级。

---

## 7. V4：结构光成像后处理

结构光只作为**探测光**，**不**加进泵浦热源。这是项目所有后处理流程的核心。

### 7.1 成像模型（标量、单波长）

probe 平面波 + 正弦条纹掩模 → 经过 plume 与坑形 → 在相机上成像：

```
I(x, y, t) = A(x, y, t) + B(x, y, t) * cos( 2π f x + phi(x, y, t) )
```

各项物理映射：

| 量 | 物理来源 | 数学映射 |
| --- | --- | --- |
| A(x, y, t) | 普通照明强度（plume 吸收 + 表面反射） | A = I0 * T_plume(x,y,t) * R_surface(x,y,t) |
| B(x, y, t) | 条纹调制度（plume 散射、表面失焦削弱） | B = A * γ(x, y, t)，γ ∈ [0, 1] |
| phi(x, y, t) | 表面位移 + plume 折射率扰动 | phi = (4π / λ) * h(x, y, t) + (2π / λ) * ∫Δn dz |

其中：

- `h(x, y, t)`：表面位移（V2 的 crater_depth 顶面投影 + V3 的 plume 反推）；
- `Δn(x, y, z, t)`：折射率扰动（由 V3 的 plume density 经 Gladstone-Dale 关系换算）。

### 7.2 输出

1. `plain_intensity_seq/`：普通照明序列（每帧一张 PNG）
2. `fringe_seq/`：结构光原始条纹序列
3. `wrapped_phase_seq/`：FFT 解调出的相位（包裹）
4. `unwrapped_phase_seq/`：解包后的相位（带 quality mask）
5. `phase_gradient_seq/`：∂phi/∂x（折射率梯度的代理）
6. `modulation_seq/`：B(x,y,t) / A(x,y,t)
7. `quality_mask_seq/`：1 = 可靠区域，0 = 被 plume 完全遮挡 / 烧蚀坑外缘
8. `phase_contrast_vs_t.csv`：定量曲线

### 7.3 关键脚本（V4 阶段实现）

- `scripts/synth_structured_light.py`：把 V2 + V3 输出合成条纹图
- `scripts/demodulate_fringe_fft.py`：FFT / Goldstein / four-step shifting 解调
- `scripts/make_four_layer_figure.py`：组装四层最终图

---

## 8. 最终四层图设计

```
+------------------------------+------------------------------+
|  Layer 1: plain intensity    |  Layer 2: structured-light  |
|  (普通照明序列)              |  fringe (条纹原始)          |
+------------------------------+------------------------------+
|  Layer 3: demodulated        |  Layer 4: quantitative       |
|  - wrapped/unwrapped phase   |  - crater depth(t)           |
|  - modulation B/A            |  - crater radius(t)          |
|  - phase gradient ∂phi/∂x    |  - plume front R(t)          |
|  - quality mask              |  - plume velocity v(t)       |
|                              |  - phase contrast / SNR      |
+------------------------------+------------------------------+
```

每个 case（不同 Ep, w₀, τ）生成一张这样的四层图，
做"参数 → 成像表现"的对比。

> 实现上：Layer 1/2 是动图（mp4 / gif），Layer 3 取关键帧，Layer 4 是曲线图。

---

## 9. 文件结构规划（目标态）

```
LS-DYNA-ablation/
├─ config/                                 # YAML 配置（V1.5 起引入）
│   ├─ experiment.yaml                     # 全局实验设定（温度阈值、单位制、约定）
│   ├─ material_silicon.yaml               # 硅物性（含 V2 起的温度依赖表）
│   ├─ laser_cases.yaml                    # V1.5 扫描矩阵（Ep, w0, tau, A, lambda）
│   └─ imaging_structured_light.yaml       # V4 成像参数（f, λ_probe, 放大率等）
│
├─ models/                                 # LS-DYNA 关键字文件
│   ├─ v1_thermal_silicon.k                # V1：纯热模型（已完成）
│   ├─ v2_threshold_ablation.k             # V2：与 V1 相同 .k，后处理在 Python
│   └─ v3_sph_ejecta.k                     # V3 路线 B 才用
│
├─ scripts/
│   ├─ _build_v1_mesh.py                   # 已有
│   ├─ run_v1.ps1                          # 已有
│   ├─ run_v1p5_sweep.ps1                  # V1.5 新增：批量跑 case
│   ├─ run_v2.ps1                          # V2 入口
│   ├─ check_v1_outputs.py                 # 已有
│   ├─ check_v1p5_outputs.py               # V1.5：批量自检
│   ├─ check_v2_outputs.py                 # V2 自检
│   ├─ plot_v1_temperature.py              # 已有：单半截面云图
│   ├─ plot_v1_full_section.py             # 已有：完整镜像截面（任务 A）
│   ├─ plot_v1_revolved_3d.py              # 已有：3D 旋转示意（任务 B）
│   ├─ extract_v2_crater.py                # V2 新增：温度阈值 → 坑形
│   ├─ synth_structured_light.py           # V4 新增：合成条纹图
│   ├─ demodulate_fringe_fft.py            # V4 新增：FFT 相位解调
│   └─ make_four_layer_figure.py           # V4 新增：组装四层图
│
├─ results/
│   ├─ v1/                                 # 已有
│   │   ├─ d3plot, tprint, messag, d3hsp, ...
│   │   └─ figures/
│   │       ├─ v1_temperature_100ns.png
│   │       ├─ v1_temperature_evolution.png
│   │       ├─ v1_temperature_full_section_100ns.png      ← 任务 A
│   │       ├─ v1_temperature_full_section_evolution.png  ← 任务 A
│   │       ├─ v1_full_section_summary.csv                ← 任务 A
│   │       └─ v1_temperature_3d_revolved_100ns.png       ← 任务 B
│   ├─ v1p5/case_XX/                       # V1.5 每个 case 独立目录
│   ├─ v2/                                 # V2 输出
│   ├─ v3/                                 # V3 plume / shock
│   └─ figures/                            # 跨版本汇总图、final 四层图
│
└─ docs/
    ├─ README_v1.md                        # 已有
    ├─ README_v2.md                        # V2 完成后补
    └─ simulation_plan_structured_light_ablation.md   # 本文件
```

---

## 10. 下一步执行建议

**坚决不要做的事**：

- ❌ 不要直接跳到 V3 SPH
- ❌ 不要直接做结构光动态仿真
- ❌ 不要在 V1 上改物理参数
- ❌ 不要在 LS-DYNA 里加 plume / 蒸气 / 相变（这些都在 Python 侧）

**推荐推进顺序**：

1. **任务 A/B 完成后**：用户人眼审阅 `results/v1/figures/` 下的镜像截面图 + 3D 旋转图，
   确认"轴对称中心热斑"的解读没有歧义；
2. **V1.5（参数扫描）**：见 §4，目的是回答"用什么 Ep + τ + λ 能真的烧硅"，
   而不是真去烧。中间产物是阈值图；
3. **V2（阈值烧蚀坑）**：见 §5，纯 Python 后处理，输出 crater(t) 曲线。
   这是给 V3 / V4 的几何输入；
4. **V3（plume / shock）**：先做路线 A（解析 + Python），再视精度需求决定要不要切到路线 B（SPH）；
5. **V4（结构光后处理）**：见 §7，把 V2 + V3 的结果合成成条纹图，做相位解调，
   产生最终四层图。

每一步完成后**必须人工审阅**：检查输出图、看 `check_*_outputs.py` 是否过、
看物理直觉是否对。**不要让 Agent 一口气把 V1.5 → V4 全做完**。

---

## 11. 下一步 Prompt 模板

把下面任一段直接复制给 Cursor Agent 即可启动下一阶段。
**注意**：每个 Prompt 都假设当前 V1 已通过、不改 V1、不改物理参数。

---

### 11.1 Prompt 1：V1.5 — 真实实验参数扫描

```
背景：
当前项目 D:\cxh-daima\LS-DYNA-ablation 已完成 V1 纯热模型（500x200 um 硅片、
高斯激光 100 ns 脉冲、峰值温度 501.5 K），见 README_v1.md。

任务：实现 V1.5 参数扫描，目的是在不引入侵蚀/相变的前提下，找出使硅表面达到
熔点 1687 K / 汽化温度 3538 K 的 (Ep, w0, tau, A, lambda) 区间。

约束：
- 不修改 V1 的 .k 文件、不修改 V1 物理参数；
- 不做 SPH、不做 plume、不做相变；
- 复用 V1 的 mesh 与材料卡，只改激光参数和（必要时）顶面网格加密；
- 所有路径用 Path(__file__).resolve().parents[1]；
- 图中文字全英文，注释用英文；
- 每个 case 独立的结果目录 results/v1p5/case_XX/；
- 单个 case 完成后 messag 必须 Normal termination 才算通过。

新增文件：
1. config/laser_cases.yaml             # 扫描矩阵（tau in [100, 30, 10, 3, 1, 0.3, 0.1] ns，
                                         w0 in [10, 20, 35, 50] um，
                                         Ep in [0.5, 5, 50, 200] uJ，
                                         A  in [0.4, 0.6]）
2. scripts/_build_v1p5_case.py         # 把每个 YAML case 编译成一份 .k
3. scripts/run_v1p5_sweep.ps1          # 顺序/并行跑全部 case
4. scripts/check_v1p5_outputs.py       # 批量自检 (Normal termination + Tmax 阈值)
5. scripts/plot_v1p5_thresholds.py     # 出 Tmax_vs_Ep / melt_map / vapor_map

输出图：
- results/v1p5/figures/Tmax_vs_Ep.png
- results/v1p5/figures/melt_threshold_map.png
- results/v1p5/figures/vapor_threshold_map.png

完成后只汇报：
- 新增了哪些文件；
- 跑了多少个 case，多少通过；
- melt / vapor 阈值在 (Ep, tau) 平面的边界在哪；
- 是否动过 V1 文件。

不要进入 V2。
```

---

### 11.2 Prompt 2：V2 — 阈值烧蚀坑后处理

```
背景：
V1.5 已经在 results/v1p5/case_XX/ 跑出多组温度场。任务是基于温度阈值
（T_melt = 1687 K, T_vap = 3538 K）在 Python 后处理里推导等效坑形，
不在 LS-DYNA 里做侵蚀 / SPH。

约束：
- 完全 Python 后处理，不重跑 LS-DYNA；
- 不改任何已有 .k 文件；
- 输入：results/v1p5/case_XX/tprint（与 V1 同格式）；
- 一个 case 一个独立子目录 results/v2/case_XX/；
- 累积型 mask：一旦某节点曾经 >= T_vap，标记为 "已剥离"，后续帧仍为 1；
- 熔融区可逆，跟随温度回落；
- 全英文图。

新增文件：
1. scripts/extract_v2_crater.py        # 单个 case：tprint → crater_mask + melt_mask
2. scripts/run_v2_batch.ps1            # 跑所有 case
3. scripts/check_v2_outputs.py         # 数值守恒检查
4. scripts/plot_v2_crater.py           # 出 crater_depth(t), crater_radius(t), 最终坑形

输出（每个 case）：
- crater_depth_t.csv
- crater_radius_t.csv
- melt_depth_t.csv
- ablation_depth_t.csv
- removed_region_map_{t}.png（关键帧）
- final_crater_section.png

汇总输出：
- results/v2/figures/crater_depth_vs_Ep.png
- results/v2/figures/crater_radius_vs_Ep.png
- results/v2/figures/ablation_volume_check.png    # 与 Ep*A/dHv/rho 对比

数值守恒检查（必须做）：
- ablated volume 与 Ep * A / dHv / rho 在 10x 量级内对齐
  （考虑到我们忽略了相变焓 ΔH_m、忽略了 plume 反推压力等，10x 是合理上限）
- 单调性：crater_depth(t)、crater_radius(t) 单调非降

完成后只汇报：
- 新增文件清单；
- 各 case 的最终坑深 / 坑径；
- 守恒检查通过率；
- 是否动过 V1 / V1.5 文件。

不要进入 V3。
```

---

### 11.3 Prompt 3：V4 — 结构光后处理 + 四层图

```
背景：
V2 已经在 results/v2/case_XX/ 给出每个 case 的 crater_depth(t), crater_radius(t),
final_crater_section.png 等。V3 还没做，先用解析模型生成 plume 等效几何
（plume_front_R(t) = R0 + v0*(t-t0)，参数从 V2 的 ablation_volume 反推）。

任务：
合成结构光条纹图，做 FFT 相位解调，输出最终四层图，逐 case 对比。

约束：
- 不重跑 LS-DYNA，不做新物理；
- 全 Python；
- 标量场，单波长（默认 800 nm），可用 imaging_structured_light.yaml 改；
- 所有图英文，300 dpi PNG；
- 一个 case 一个 results/v4/case_XX/ 子目录。

物理模型：
I(x,y,t) = A(x,y,t) + B(x,y,t) * cos(2*pi*f*x + phi(x,y,t))
其中：
- A = I0 * T_plume * R_surface
- B = A * gamma_plume      (gamma in [0, 1], plume 散射造成调制度衰减)
- phi = (4*pi/lambda) * h(x,y,t) + (2*pi/lambda) * integral(delta_n) dz

输入：
- V2: crater_depth(t), crater_radius(t), final_crater_section
- V3-light（解析）: plume_front_R(t), plume_density_proxy_map

新增 config:
- config/imaging_structured_light.yaml
    fringe_freq_lp_per_mm, lambda_probe_nm, n_frames, field_of_view_um,
    camera_pixel_um, magnification, gladstone_dale_coefficient

新增脚本：
1. scripts/synth_plume_analytic.py     # 从 V2 ablation_volume → V3 plume 解析参数
2. scripts/synth_structured_light.py   # 合成 I(x,y,t)；同时输出 plain 和 fringe 序列
3. scripts/demodulate_fringe_fft.py    # FFT 解调 wrapped_phase / modulation / quality
4. scripts/make_four_layer_figure.py   # 拼装四层图

输出（每个 case）：
- plain_intensity_seq/         frame_NNN.png
- fringe_seq/                  frame_NNN.png
- wrapped_phase_seq/           frame_NNN.png
- unwrapped_phase_seq/         frame_NNN.png
- phase_gradient_seq/          frame_NNN.png
- modulation_seq/              frame_NNN.png
- quality_mask_seq/            frame_NNN.png
- phase_contrast_vs_t.csv
- four_layer_figure.png        ← 最终汇报图（4 块拼接）

约束检查：
- wrapped_phase 在 plume 遮蔽区应自动被 quality_mask 标记；
- crater_depth 与 解调出的 h(x=0, y=0, t) 应该在 50% 误差内对齐；
- modulation 在 plume 浓密区应 < 0.2；
- 四层图必须包含定量曲线 (crater_depth, crater_radius, plume_front_R, plume_v, phase_contrast)。

完成后只汇报：
- 新增文件清单；
- 各 case 的关键定量数：max crater_depth、max plume_front、min modulation；
- 是否动过 V1/V1.5/V2 文件。
```

---

## 12. 附录：参数 / 常数速查

| 量 | 符号 | 数值 | 单位 |
| --- | --- | --- | --- |
| 硅密度 | ρ | 2330 | kg/m³ |
| 硅比热 | cp | 700 | J/(kg·K) |
| 硅导热 | k | 150 | W/(m·K) |
| 硅热扩散率 | α = k/(ρ cp) | 9.2×10⁻⁵ | m²/s |
| 硅熔点 | T_melt | 1687 | K |
| 硅沸点 | T_vap | 3538 | K |
| 硅蒸发焓 | ΔH_v | 1.79×10⁷ | J/kg |
| 真空光速 | c | 2.998×10⁸ | m/s |
| Gladstone-Dale (Si plume) | K_GD | ~1.5×10⁻⁴ | m³/kg（待校准） |

单位换算：

- 1 J/cm² = 10⁴ J/m²
- 1 mJ on 70 µm 直径光斑 ≈ 26 J/cm²
- LS-DYNA mm-ms-kg-kN-GPa 制：1 kW/mm² = 10⁹ W/m²

---

*文档结束。本方案在 V1.5 / V2 完成后会回填实测数据并修订。*
