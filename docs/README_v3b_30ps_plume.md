# V3B — 30 ps 参数化羽流 / 冲击波 / 抛射物模型

V3B 在 **V2.6 温度阈值等效 crater driver**（`ep_5p0uj`）基础上，建立**分析型 / 参数化** plume、shock front 与 ejecta 展示模型，为后续 **V4 结构光** 与 **LS-PrePost 可视化** 路线做准备。

**不是** LS-DYNA 显式多相喷射，**不是** SPH / ALE，**不修改** V1.7 / V2.6 原始结果。

---

## 1. 目的

- 将 V2.6 等效烧蚀坑几何与时刻，转译为可扩展的 plume / shock 参数序列；
- 生成 density proxy、plain imaging、ejecta 粒子 proxy 等**展示型**图；
- 输出机器可读 CSV / JSON，供 V4 结构光条纹合成读取。

---

## 2. 为什么使用 V2.6 driver

V2.6 从 V1.7 30 ps 局部细网格 `tprint` 提取 threshold-equivalent melt / vapor / crater 指标。  
主 case **`ep_5p0uj`** 是唯一 `vapor/ablation candidate`，作为 V3B 驱动。

Driver 路径：

```
results/v26_30ps_threshold_ablation/v3b_driver/v3b_driver_ep5p0uj.json
```

由 `scripts/export_v3b_driver_from_v26.py` 导出（只读 V2.6 metrics）。

---

## 3. t0 为何取 max crater-depth frame（50.929 ps）

| 参考帧 | t (ps) | 物理含义 |
| --- | ---: | --- |
| Tmax peak | 30.03 | 最高 lattice 温度、较大坑半径 |
| **max crater depth** | **50.93** | 全时间最大等效坑深 |

V3B **plume launch time** `t0_ps` 默认使用 driver 的 `recommended_t0_ps = 50.929 ps`（max crater-depth frame），因为**几何尺度**（坑深 / 坑半径）在该时刻达到 V2.6 定义的全局最大坑深，更适合作为羽流“起爆”几何输入。

Tmax peak frame（30.03 ps）仍保留在 driver 中，用于**热峰值 / 能量说明**，不作为默认 plume 几何 t0。

---

## 4. 模型公式（t_after = max((t_ps - t0_ps)/1000, 0) ns）

**Plume front**

```
R_plume = R0 + v0 * t_after^b
```

**Shock front (V3B.1)**

```
R_shock_raw = R0 + c0 * dt + beta * sqrt(dt)
R_shock = max(R_shock_raw, R_plume + shock_ahead_min)
```

**Density proxy (V3B.1 semi-ellipsoid + root term)**

```
n0 = 1 / (1 + dt / density_decay_ns)^p
sigma_r = max(sigma_r_min, frac_r * R_plume)
sigma_z = max(sigma_z_min, frac_z * R_plume)
z_center = z_center_frac * R_plume
n(r,z) = n0 * exp(-(r/sigma_r)^2 - ((z-z_center)/sigma_z)^2) + 0.25*n0*root(r,z)   (z >= 0)
root = exp(-(r/(0.8*R0))^2) * exp(-(z/10)^2)
```

**Ejecta proxy**

```
ejecta_height = h_max * (1 - exp(-t_after/0.2)) * exp(-t_after/decay_ns)
ejecta_radius = R0 + radial_spread * (1 - exp(-t_after/1.0))
```

**坐标系（展示用）**  
表面 `z = 0`，向上为 plume / shock 区域；与 V1.7/V2.6 的 `z = 5 µm` 硅片网格坐标**独立**。

---

## 5. 文件与脚本

| 路径 | 说明 |
| --- | --- |
| `config/plume_model_v3b_30ps.toml` | 模型参数 |
| `scripts/build_v3b_plume_model.py` | 生成 metrics CSV |
| `scripts/plot_v3b_plume_shock.py` | 出图 |
| `results/v3b_30ps_plume/v3b_plume_shock_metrics.csv` | 0–5 ns 时间序列 |
| `results/v3b_30ps_plume/v3b_driver_used.json` | 本次运行使用的 driver 快照 |
| `results/v3b_30ps_plume/v3b_t0.txt` | t0_ps / t0_ns |
| `results/v3b_30ps_plume/figures/` | PNG 图 |

---

## 6. 运行

```powershell
.\.venv\Scripts\python.exe scripts\export_v3b_driver_from_v26.py
.\.venv\Scripts\python.exe scripts\build_v3b_plume_model.py
.\.venv\Scripts\python.exe scripts\plot_v3b_plume_shock.py
```

---

## 7. 哪些图可进入 V4 结构光

| 图 | V4 用途 |
| --- | --- |
| `v3b_density_proxy_sequence.png` | 条纹调制度衰减 / 相位畸变输入场 |
| `v3b_plain_imaging_sequence.png` | plain-light 对照基线 |
| `v3b_plume_shock_front_vs_time.png` | 前沿时刻标定 |
| `v3b_ejecta_particle_proxy_sequence.png` | 散射体分布示意（非物理粒子） |
| `v3b_driver_summary.png` | 汇报 / PPT 总览 |
| `v3b_early_time_zoom_sequence.png` | 0–500 ps 早期 plume/ejecta 细节 |
| `v3b_v4_input_fields.png` | V4 输入场检查（density / B / phi proxy） |

V4 应读取 `v3b_plume_shock_metrics.csv` + `v3b_driver_used.json` + `density_field()`（或等价函数），而非硬编码常数。

---

## 8. 局限性

- 非真实多相 LS-DYNA；无单元删除 / 蒸发质量守恒；
- 非 SPH / ALE；ejecta 为固定 seed 的**展示型**粒子 proxy；
- plume density 为 Gaussian proxy，非 CFD / 辐射传输；
- R0、d0 来自 **threshold-equivalent** 定义（T ≥ T_vap），非实测 crater 形貌。

---

## 9. 后续升级路线

| 阶段 | 内容 |
| --- | --- |
| **V3C** | ejecta / plume proxy → LS-PrePost 可读展示格式 |
| **V4** | 结构光条纹调制度衰减、相位畸变合成 |
| 后期 | 真实 SPH/ALE 或 LS-DYNA 单元侵蚀（超出当前 V3B 范围） |

---

## 10. V2.6 → V3B driver handoff（摘要）

1. **Tmax peak frame** 与 **max crater-depth frame** 时刻不同、物理含义不同。  
2. ep_5p0uj：**30.03 ps** → Tmax ≈ 5604 K，crater d ≈ 0.031 µm，r ≈ 17.39 µm；**50.93 ps** → 最大坑深 d ≈ 0.042 µm，r ≈ 15.03 µm。  
3. V3B 几何 t0 用 **50.93 ps**；热峰值说明用 **30.03 ps**。  
4. V2.6 / V3B 均为 **threshold-equivalent postprocess**，不是显式多相 / 喷射 LS-DYNA 求解。

---

## 11. V3B.1 visual calibration

### 为什么要调整

原 V3B 参数下，0–5 ns 内 `R_plume` 仅约 40 µm、`R_shock` 后期被 clamp 后与 plume 几乎重合，density proxy 贴近表面、plain imaging 暗柱不明显。  
**V3B.1 只做展示型参数校准**，使羽流 / 冲击波 / ejecta 在侧视图中更易辨认，**不改变** V1.7 / V2.6 真实热结果或 driver 数值。

### 物理定位（仍为 proxy）

| 量 | 定位 |
| --- | --- |
| density | plume density **proxy**（半椭球 Gaussian + 表面 root 项） |
| ejecta particles | **display-only** proxy，非质量守恒 SPH 粒子 |
| shock front | 解析前沿 `R_shock(t)`，非 CFD 激波 |
| plain imaging | Beer–Lambert 风格 extinction 合成 |

不能把粒子点云当作真实喷出质量；不能把 density proxy 当作 CFD 数密度。

### V3B.1 关键参数（`config/plume_model_v3b_30ps.toml`）

| 参数 | 值 |
| --- | ---: |
| `v0_um_per_ns_power` | 35.0 |
| `b` | 0.72 |
| `c0_um_per_ns` | 18.0 |
| `beta_um_per_sqrt_ns` | 25.0 |
| `shock_ahead_min_um` | 5.0 |
| `density_decay_ns` | 2.0 |
| `density_decay_power` | 1.2 |
| `tau_eff` (plain imaging) | 2.5 |
| `max_ejecta_height_um` | 120.0 |
| `radial_spread_um` | 180.0 |
| `decay_ns` (ejecta) | 2.5 |

目标：5 ns 时 `R_plume` ≈ 120–140 µm，`R_shock` ≈ 180–220 µm，shock 始终领先 plume ≥ `shock_ahead_min_um`。

### V4 将使用

- `results/v3b_30ps_plume/v3b_plume_shock_metrics.csv`
- `results/v3b_30ps_plume/v3b_driver_used.json`
- `results/v3b_30ps_plume/figures/v3b_v4_input_fields.png`（输入场检查）
- `scripts/plot_v3b_plume_shock.py` 内 `density_field()`（或 build 侧等价逻辑）

### 后期 LS-PrePost / 真实喷射路线

| 阶段 | 内容 |
| --- | --- |
| **V5A** | 温度阈值单元删除 / erosion → 烧蚀坑 d3plot（LS-PrePost 可看） |
| **V5B** | SPH/ALE 或等效喷射粒子，展示材料喷出与粒子云 |
| **V4** | 结构光条纹调制度衰减、相位畸变（基于 V3B.1 proxy 场） |
