# V1.7 — Local refined 2D axisymmetric 30 ps thermal model

V1.7 builds a **local, non-uniform, fine-surface-mesh** LS-DYNA
single-temperature thermal model for 30 ps picosecond laser heating on
silicon.  It is the natural successor to V1.6's scale analysis.

V1.7 does **not** modify V1 / V1.5 / V2 / V3A / V1.6 files or results.
It does **not** run element erosion, SPH, plume, structured light, or
real 3-D models.

> **Important:** V1.7 is a *single-temperature lattice* LS-DYNA
> approximation.  Compare its temperatures with the **lattice branch
> `Tl`** of the V1.6 TTM prototype — **not** with `Te` and not as a
> complete ultrafast ablation model.

Later stages may add LS-DYNA element deletion and SPH/ALE ejecta
visualisation routes; **V1.7 does not include those**.

---

## 1. Why a local refined mesh?

V1.6 showed that at 30 ps the thermal diffusion length is only
**L_diff ≈ 53 nm**, while V1 uses **dz = 2500 nm** at the surface
(~**235× too coarse**).  A full wafer (r = 500 µm, z = 200 µm) with
10 nm surface cells would require billions of elements.

V1.7 therefore uses a **local domain**:

| quantity | V1 / V1.5 | V1.7 local |
| --- | --- | --- |
| radius | 500 µm | **100 µm** |
| thickness | 200 µm | **5 µm** |
| surface dz | 2500 nm | **10 nm** |
| pulse | 100 ns | **30 ps** |

This resolves the heat-affected zone while keeping the model runnable
on a workstation.

---

## 2. Mesh specification (built mesh)

From `models/v17_30ps_local/v17_mesh_summary.csv`:

| parameter | value |
| --- | --- |
| r range | 0 – **100 µm** |
| z range | 0 – **5 µm** (top surface = laser side) |
| NR × NZ | **80 × 92** |
| nodes | **7 533** |
| shell elements | **7 360** |
| dr | **0.5 – 2.5 µm** (fine near axis) |
| dz (min at top surface) | **10.0 nm** |
| dz (max in bulk) | **465 nm** |
| axis-top node id | **7453** (r = 0, z = top) |

Z-grid strategy (bottom → top):

- Bottom bulk (z = 0 → ~3 µm): coarsening steps 100 → 465 nm
- Mid layer (~3 → 4.5 µm): ~50–100 nm steps
- Top **0.5 µm** below surface: **10 nm** constant (matches V1.6 L/5 ≈ 10.6 nm)

R-grid strategy:

- r = 0 – 10 µm: dr = **0.5 µm**
- r = 10 – 50 µm: dr = **1.0 µm**
- r = 50 – 100 µm: dr = **2.5 µm**

Node coordinates are exported to `v17_mesh_nodes.csv` for post-processing.

---

## 3. Why still 2D axisymmetric?

The experiment uses a circular Gaussian spot on a flat wafer — good
rotational symmetry.  2D axisymmetric ELFORM=15 shells give the same
physics as a 3-D disk at **~1/1000th** the element count.

As explained in V1.6 (`v16_axisymmetric_coordinate_explanation_zoomed.png`),
the hot zone appears in the **upper-left corner** of the raw LS-PrePost
view because r = 0 is the left edge.  Mirror post-processing places it
at the **top centre** of the full cross-section.

---

## 4. Why not 3-D?

At 30 ps with 10 nm surface cells, a full 3-D disk (r = 500 µm,
z = 200 µm) would need **10⁹–10¹²** elements — impractical on a student
workstation.  See `docs/README_3d_feasibility_after_v17.md` (planned)
for the three optional 3-D routes.  V1.7 intentionally stays 2-D.

---

## 5. Why V1.7 is not full TTM

LS-DYNA's built-in thermal solver is **single-temperature**.  At 30 ps
the electron and lattice subsystems are not in equilibrium; electrons
heat first and transfer energy to the lattice over τ_ep ≈ 10–100 ps.

V1.6's 1D TTM prototype (`scripts/ttm_1d_30ps_prototype.py`) models
this explicitly.  V1.7 deposits the laser energy directly into the
**lattice** with a Gaussian temporal pulse — it will generally
**over-estimate** lattice heating rate compared to TTM at the same Ep,
or require higher Ep to reach the same Tl.

---

## 6. Comparison with V1.6 TTM (lattice Tl)

| Ep (µJ) | V1.6 TTM Tl_peak (prototype) | expected V1.7 behaviour |
| ---: | ---: | --- |
| 0.5 | 828 K | below melt |
| 1.0 | 1 348 K | near melt |
| 2.0 | 2 364 K | above melt |
| 5.0 | 5 237 K | above vapor |

Single-temperature analytic estimates (V1.6) give Ep_melt ≈ 0.42 µJ —
much lower than TTM because they ignore electron–phonon lag.

After running V1.7, compare `v17_case_summary.csv` T_peak against the
TTM table above.  Agreement within ~2× is acceptable given the
single-temperature approximation; larger gaps indicate the need for
coupled TTM–FD or shorter ITS tuning.

---

## 7. Selected cases

| case | Ep (µJ) | I0 peak (kW/mm²) | V1.6 TTM note |
| --- | ---: | ---: | --- |
| `ep_0p5uj` | 0.5 | 1.73×10³ | below TTM melt |
| `ep_1p0uj` | 1.0 | 3.46×10³ | near TTM melt |
| `ep_2p0uj` | 2.0 | 6.91×10³ | above TTM melt |
| `ep_5p0uj` | 5.0 | 1.73×10⁴ | above TTM vapor |

All cases: τ = 30 ps, w0 = 35 µm, A = 0.5, t_end = 5 ns.

Output snapshot targets (for plotting): 0, 10, 30, 60, 100, 300 ps,
1 ns, 2 ns, 5 ns.

### LS-DYNA database output interval

The `[output_times_ns]` list in the TOML is for **post-processing targets only**.
The actual tprint / d3plot times are controlled by the constant **DT** on these
keyword cards (all set to the same value in `_build_v17_case.py`):

| keyword | role |
| --- | --- |
| `*DATABASE_BINARY_D3PLOT` | binary state dumps |
| `*DATABASE_BINARY_D3THDT` | binary thermal state |
| `*DATABASE_TPRINT` | ASCII node temperatures (used by check/plot scripts) |
| `*DATABASE_GLSTAT` | global energy balance |

**Previous (bug):** `DT = t_end / n_plots_desired = 5 ns / 40 = 0.125 ns`
(~125 ps).  Snapshots landed at 0, 0.13, 0.25, … ns — missing 10/30/60/100 ps.

**Current:** `[database_output] dt_output_ns = 0.01` → `DT = 1.0×10⁻⁸ ms`
(0.01 ns = 10 ps).  ~500 frames over 5 ns; hits all early-pulse targets.

LS-DYNA does not support piecewise DT on these cards in our setup, so a
uniform 0.01 ns interval is used (finest segment of the reference schedule
0–0.2 ns / 0.2–1 ns / 1–5 ns documented in the TOML).

To validate output timing without a full 5 ns rerun, rebuild and run only
`ep_1p0uj` — existing results under `results/v17_30ps_local/ep_1p0uj/` are
left untouched until you explicitly rerun.

---

## 8. Files

| path | role |
| --- | --- |
| `config/v17_30ps_local_mesh.toml` | all parameters |
| `scripts/_build_v17_case.py` | internal mesh + .k builder |
| `scripts/build_v17_30ps_local_mesh.py` | generate mesh + 4 .k files |
| `scripts/run_v17_selected.ps1` | LS-DYNA launcher |
| `scripts/check_v17_outputs.py` | automated checks |
| `scripts/plot_v17_30ps_temperature.py` | 6 PNG figures |
| `models/v17_30ps_local/v17_*.k` | LS-DYNA inputs |
| `models/v17_30ps_local/v17_case_registry.csv` | case registry |
| `models/v17_30ps_local/v17_mesh_summary.csv` | mesh statistics |
| `models/v17_30ps_local/v17_mesh_nodes.csv` | node coordinates |
| `results/v17_30ps_local/<case>/` | LS-DYNA outputs |
| `results/v17_30ps_local/v17_case_summary.csv` | check summary |
| `results/v17_30ps_local/figures/` | post-process PNGs |

---

## 9. How to run

```powershell
# Step 1 -- build mesh + .k (re-run after changing database_output in TOML)
.\.venv\Scripts\python.exe scripts\build_v17_30ps_local_mesh.py

# Step 2 -- preview what will run
.\scripts\run_v17_selected.ps1 -DryRun

# Step 3 -- run ONE case first (recommended; 0.01 ns output => ~500 frames / case)
.\scripts\run_v17_selected.ps1 -OnlyCase ep_1p0uj -Ncpu 4

# Step 4 -- run all selected cases
.\scripts\run_v17_selected.ps1 -Ncpu 4

# Step 5 -- check + plot
.\.venv\Scripts\python.exe scripts\check_v17_outputs.py
.\.venv\Scripts\python.exe scripts\plot_v17_30ps_temperature.py
```

Each .k file is ~816 KB ASCII.  With 7533 nodes the thermal solve is
much lighter per step than V1, but the explicit timestep (ITS ≈ 0.05 τ
≈ 1.5 ps) over 5 ns implies many steps — expect **minutes to hours**
per case depending on hardware.

---

## 10. When to enter V2.6

Proceed to **V2.6** (30 ps threshold ablation post-processing) when:

1. At least **`ep_2p0uj`** and **`ep_5p0uj`** show `Normal termination`.
2. `check_v17_outputs.py` reports `peak_at_axis_top = yes`.
3. T_peak ordering is monotonic in Ep (0.5 < 1.0 < 2.0 < 5.0 µJ).
4. Qualitative agreement with V1.6 TTM Tl trends (melt near 1–2 µJ,
   vapor near 2–5 µJ).

V2.6 will mirror V2's threshold metrics but read V1.7 tprint + the
non-uniform `v17_mesh_nodes.csv` grid.

---

## 11. What V1.7 does NOT change

- `models/v1_thermal_silicon.k`, `models/v15_cases/**` — unchanged
- `results/v1/**`, `results/v15/**`, `results/v2/**`, `results/v3a/**`,
  `results/v16_30ps/**` — unchanged
- `config/laser_cases_v15.toml`, `config/plume_model_v3a.toml`,
  `config/v16_30ps.toml` — unchanged

All V1.7 writes go to `models/v17_30ps_local/` and `results/v17_30ps_local/`.
