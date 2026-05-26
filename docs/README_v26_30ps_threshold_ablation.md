# V2.6 — 30 ps 温度阈值等效烧蚀坑后处理

V2.6 在 **V1.7 局部细网格 30 ps LS-DYNA 温度场** 基础上，用 Python 做**温度阈值等效烧蚀坑**分析。  
**不是** LS-DYNA 单元删除、**不是** SPH/ALE、**不是** 真实 3D 热-力耦合。

> 驱动模型：V1.7 local mesh  
> 主 case：`ep_5p0uj`（唯一 `vapor/ablation candidate`）

---

## 1. 物理意义

| 区域 | 条件 | 含义 |
| --- | --- | --- |
| melt zone | T ≥ T_melt (1687 K) | 等效熔化区 |
| vapor / ablation zone | T ≥ T_vap (3538 K) | 等效蒸气/烧蚀区 |
| crater | depth = vapor_depth, radius = vapor_radius | 等效坑几何 |

初始温度 T_init = 300 K。  
深度：从上表面 (z = 5 µm) 沿轴线 r = 0 向下，线性插值求阈值交点。  
半径：在上表面沿 r 向外，线性插值求阈值交点。

V1.7 为 2D 轴对称 half-section；LS-PrePost 默认视图中 r=0 在**左侧**、激光点在**左上角**。  
V2.6 出图时对 r 做镜像，使热区/烧蚀坑显示在 **top centre**。

---

## 2. 输入

| 路径 | 说明 |
| --- | --- |
| `models/v17_30ps_local/v17_case_registry.csv` | case 列表 |
| `models/v17_30ps_local/v17_mesh_nodes.csv` | 非均匀 r-z 节点坐标 |
| `results/v17_30ps_local/<case>/tprint` | V1.7 温度输出 |
| `config/v17_30ps_local_mesh.toml` | T_melt / T_vap |

已运行 case（2026-05）：`ep_1p0uj`, `ep_2p0uj`, `ep_5p0uj`；`ep_0p5uj` 可缺失。

---

## 3. 输出

```
results/v26_30ps_threshold_ablation/
  <case>/v26_threshold_metrics.csv    # 每个时间步指标
  v26_case_summary.csv                # case 汇总
  figures/
    v26_case_comparison_table.png
    v26_crater_depth_vs_time.png
    v26_crater_radius_vs_time.png
    v26_crater_depth_vs_time_early.png      # 0--120 ps zoom
    v26_crater_radius_vs_time_early.png
    v26_ep5p0uj_v3b_driver_summary.png      # 2x2 V3B handoff figure
    v26_final_threshold_maps.png
    v26_ep5uj_crater_evolution.png
    v26_ep5uj_crater_evolution_zoomz.png   # surface-near z zoom (4.94--5.005 um)
    v26_ep5uj_peak_crater_profile_zoom.png   # Tmax peak frame
    v26_ep5uj_max_crater_profile_zoom.png    # max crater-depth frame
    v26_ep5uj_revolved_crater_schematic.png
```

### 指标 CSV 列

`time_ns`, `time_ps`, `Tmax_K`, `melt_exists`, `vapor_exists`,  
`melt_depth_um`, `melt_radius_um`, `vapor_depth_um`, `vapor_radius_um`,  
`crater_depth_um`, `crater_radius_um`, `r_at_Tmax_um`, `z_at_Tmax_um`

### final_regime（case 级）

- Tmax < T_melt → `no melt`
- T_melt ≤ Tmax < T_vap → `melt only`
- Tmax ≥ T_vap → `vapor/ablation candidate`

---

## 4. 脚本

| 脚本 | 作用 |
| --- | --- |
| `scripts/extract_v26_30ps_threshold_ablation.py` | 读 tprint，写 CSV |
| `scripts/plot_v26_30ps_threshold_ablation.py` | 读 CSV + tprint，出图 |
| `scripts/export_v3b_driver_from_v26.py` | 导出 ep_5p0uj V3B driver package |

复用：`parse_tprint`（V1）、`load_mesh_nodes` / `tprint_to_grid`（V1.7 check）。

---

## 5. 运行

```powershell
.\.venv\Scripts\python.exe scripts\extract_v26_30ps_threshold_ablation.py
.\.venv\Scripts\python.exe scripts\plot_v26_30ps_threshold_ablation.py
.\.venv\Scripts\python.exe scripts\export_v3b_driver_from_v26.py
.\.venv\Scripts\python.exe scripts\check_v17_outputs.py
```

**无需重跑 LS-DYNA。**

---

## 6. V2.6 → V3B driver handoff

### 两个参考帧，物理含义不同

| 参考帧 | 时刻 (ps) | Tmax (K) | crater d (µm) | crater r (µm) | 用途 |
| --- | ---: | ---: | ---: | ---: | --- |
| **Tmax peak frame** | 30.03 | 5604 | 0.031 | 17.39 | 最高温度、较大坑半径、热驱动说明 |
| **max crater-depth frame** | 50.93 | 4982 | **0.042** | 15.03 | 全时间最大等效坑深、几何尺度 |

二者**不是同一时刻**：30.03 ps 温度最高但坑尚较浅；50.93 ps 坑深最大但温度已回落。

### V3B 建议使用方式

| 用途 | 推荐参考 |
| --- | --- |
| 几何尺度（坑深/坑半径） | `max_crater_depth_frame` → t₀ ≈ 50.93 ps, d ≈ 0.042 µm, r ≈ 15.03 µm |
| 热驱动 / 能量峰值说明 | `tmax_peak_frame` → 30.03 ps, Tmax ≈ 5604 K |

机器可读 package：

```
results/v26_30ps_threshold_ablation/v3b_driver/
  v3b_driver_ep5p0uj.json
  v3b_driver_ep5p0uj.csv
  v3b_driver_ep5p0uj_readme.md
```

由 `scripts/export_v3b_driver_from_v26.py` 从 V2.6 metrics **只读**导出，不修改 CSV 数值。

### 重要限制

V2.6 仍是 **temperature-threshold equivalent postprocess**（T ≥ T_melt / T_vap 等效区），**不是**显式多相、单元删除或喷射 SPH/ALE 求解。V3B 若做羽流/冲击 analytic 模型，应将此 package 视为几何与时刻的**驱动输入**，而非 LS-DYNA 真实 ablation 边界。

---

## 7. 与 V1.7 / V3B / LS-PrePost 的关系

| 阶段 | 关系 |
| --- | --- |
| V1.7 | 提供 tprint 与非均匀网格；V2.6 只读 |
| V2 (100 ns) | 同一套阈值逻辑，V2.6 适配 30 ps + 非均匀网格 |
| V3B（计划） | 可用 V2.6 `ep_5p0uj` crater 指标作分析/羽流驱动 |
| LS-PrePost | V2.6 镜像图用于论文；LS-PrePost 仍看 half-section 左上角 |

---

## 8. 不改动的内容

- V1 / V1.5 / V1.6 / V1.7 的 `.k`、物理参数、已有 results
- 无 SPH、无真实 3D LS-DYNA、无结构光

所有 V2.6 写入 `results/v26_30ps_threshold_ablation/` 及 `docs/README_v26_30ps_threshold_ablation.md`。
