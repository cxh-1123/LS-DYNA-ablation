# V1.6 — 30 ps analytic scale + 1D TTM prototype

V1.6 is the **pre-flight stage** for the 30 ps picosecond branch of
the project.  It is **pure Python** and runs in seconds.  It does
**not** run LS-DYNA, it does **not** modify any V1 / V1.5 / V2 / V3A
file, and it does **not** introduce SPH, element erosion, plume CFD,
or structured-light demodulation.

V1.6 answers seven questions before any 30 ps LS-DYNA model is built:

1. What is the thermal diffusion length at 30 ps, 100 ps, 1 ns, 100 ns?
2. Is the current V1 mesh (`dz = 2.5 µm`) fine enough for 30 ps?
3. What surface `dz` is recommended?
4. Roughly what pulse energy is needed to reach melt / vapor at 30 ps?
5. Does a 30 ps shot need a two-temperature model (TTM)?
6. In the V1.5 TTM prototype, do `Te` and/or `Tl` reach the thresholds?
7. Why does the heated zone appear in the **upper-left corner** in
   LS-PrePost, and where does it really sit in the wafer?

---

## 1. Why not just reuse V1 with a 30 ps pulse?

V1 was built around `dz = 2.5 µm = 2500 nm` at the surface.  The
thermal diffusion length at 30 ps in silicon is

```
L_diff(30 ps) = sqrt(alpha * tau) = sqrt(9.4e-5 m^2/s * 30e-12 s) ≈ 53 nm
```

so `V1 dz / L_diff = 47` and `V1 dz / (L_diff / 5) = 235`.  The V1
mesh is two and a half orders of magnitude too coarse to resolve the
30 ps heat-affected zone.  Even at 1 ns the V1 mesh is still ~41×
too coarse, and only at 100 ns does it become ~4× too coarse (which
is why V1 / V1.5 already saw small sub-cell artefacts).

The other reason is physics: at 30 ps the electron and lattice
subsystems are not in equilibrium.  A single-temperature LS-DYNA model
overestimates the lattice temperature because it ignores the
electron–phonon coupling delay.  See §5 below.

---

## 2. Physical differences between 30 ps and 100 ns

| Aspect                          | 100 ns (V1.5)                 | 30 ps (V1.6)                                  |
| ------------------------------- | ----------------------------- | --------------------------------------------- |
| L_diff                          | ~3 µm                         | ~53 nm                                        |
| Recommended dz                  | ~600 nm                       | ~10 nm                                        |
| V1 dz / recommended dz          | ~4× too coarse                | ~235× too coarse                              |
| `Ep_melt` (1D half-inf, w0=35 µm, A=0.5) | ~24 µJ              | ~0.42 µJ                                      |
| `Ep_vap` (1D half-inf)          | ~56 µJ                        | ~0.97 µJ                                      |
| Electron–lattice equilibration  | irrelevant (tau >> tau_ep)    | dominant (tau ≈ tau_ep, TTM required)         |
| Phase-explosion / spallation    | secondary                     | possibly primary                              |
| Plume formation timescale       | ns                            | ps                                            |

The Ep estimates use a one-temperature 1D half-infinite Stefan-like
formula

```
dT_peak = 2 F_abs / sqrt(pi tau rho cp k)
F_abs   = A * 2 Ep / (pi w0^2)
=> Ep   = dT * pi * w0^2 * sqrt(pi tau rho cp k) / (4 A)
```

For 100 ns it predicts Ep_melt ≈ 24 µJ, which matches V1.5's measured
melt onset Ep ≈ 25 µJ within 5 %.  This calibrates the formula well
enough for use as a 30 ps starting estimate.

---

## 3. Diffusion length and recommended `dz`

Exact analytic values used by `estimate_v16_30ps_scales.py`:

| tau     | L_diff   | recommended dz (= L/5) | V1 dz / rec dz | adequate? | TTM regime                     |
| ------: | -------: | ---------------------: | -------------: | :-------: | ------------------------------ |
|   30 ps |   53.1 nm|              10.6 nm   |        235×    |     no    | required (tau < 100 ps)        |
|  100 ps |   97.0 nm|              19.4 nm   |        129×    |     no    | marginal (100 ps ≤ tau < 1 ns) |
|    1 ns |    307 nm|              61.3 nm   |         41×    |     no    | not required (tau ≥ 1 ns)      |
|  100 ns |   3.07 µm|               613 nm   |          4×    |     no    | not required                   |

See `results/v16_30ps/v16_30ps_scale_summary.csv` for the machine-
readable form, and `results/v16_30ps/figures/v16_diffusion_length_vs_tau.png`
and `v16_mesh_resolution_warning.png` for the plots.

> Take-away: for V1.7 the surface dz must be **5 – 20 nm**.  Below
> 5 nm we'd be over-resolving heat at the cost of LS-DYNA time-stepping;
> above 20 nm we are starting to miss the 30 ps temperature gradient.

---

## 4. 30 ps `Ep_melt` and `Ep_vap` analytic estimates

For `w0 = 35 µm`, `A = 0.5`:

| tau     | Ep_melt (1D analytic) | Ep_vap (1D analytic) |
| ------: | --------------------: | -------------------: |
|   30 ps |             0.42 µJ  |             0.97 µJ  |
|  100 ps |             0.76 µJ  |             1.78 µJ  |
|    1 ns |             2.41 µJ  |             5.62 µJ  |
|  100 ns |             24.1 µJ  |             56.2 µJ  |

These are **single-temperature** estimates.  In real picosecond
ablation the lattice temperature lags the electron temperature, so
the actual `Ep_melt` and `Ep_vap` for the *lattice* are typically a
few times higher (see §6).  The single-T estimates are still useful
as the lower bound and as a sanity-check against the V1.7 LS-DYNA
results.

Plot: `results/v16_30ps/figures/v16_threshold_energy_vs_tau.png`.

---

## 5. Why a two-temperature model is needed at 30 ps

At nanosecond timescales the electron–phonon equilibration time
`tau_ep ≈ 10 – 100 ps` for silicon is much shorter than the pulse,
so the single-temperature equation `T_lattice = T_electron = T(z, t)`
is correct.  At 30 ps these two timescales are *comparable*: the
electrons absorb essentially all of the laser energy first (their
heat capacity is small), then transfer it to the lattice over `tau_ep`.

Standard two-temperature equations (Anisimov-Kapeliovich-Perelman):

```
Ce(Te) ∂Te/∂t = ∂/∂z [ ke ∂Te/∂z ] - G (Te - Tl) + S(z, t)
Cl     ∂Tl/∂t = ∂/∂z [ kl ∂Tl/∂z ] + G (Te - Tl)
S(z, t)       = A * q_surface(t) / delta_abs * exp(-z / delta_abs)
q_surface(t)  = F_peak / (sqrt(2 pi) tau) * exp(-(t - t_c)^2 / (2 tau^2))
```

In V1.6 the equations are integrated on a 1D `z` grid (`z = 0 .. 500 nm`,
`dz = 5 nm`) with forward Euler and a CFL-throttled `dt`.

### 5.1 Prototype parameter caveats

| parameter             | value used      | reality (Si)                                          |
| --------------------- | --------------- | ----------------------------------------------------- |
| `Ce(Te) = gamma_e Te` | gamma_e = 100   | Si has a band gap; Ce is non-linear in Te             |
| `ke`                  | 150 W/(m·K)     | strongly Te-dependent; depends on free-carrier density|
| `G`                   | 1e17 W/(m³·K)   | reported values 1e16 – 5e18 W/(m³·K)                  |
| `delta_abs`           | 5 / 10 / 1000 nm| depends on doping, temperature, free-carrier abs.     |

These are **placeholders**.  The V1.6 TTM is therefore a
*qualitative scale-finder* — it shows the right behaviour (Te peaks
during the pulse, Tl follows with a delay, peak Tl is below peak Te,
the equilibration time-scale is ~ 50 – 100 ps) — but its Ep_melt /
Ep_vap predictions are not calibrated.  All plots and CSVs from
`ttm_1d_30ps_prototype.py` carry the `PROTOTYPE PARAMETERS` flag.

### 5.2 Prototype results (V1.6 default: tau=30 ps, lambda=266 nm)

| Ep (µJ) | Te_peak (K) | Tl_peak (K) | lattice melt? | lattice vapor? |
| ------: | ----------: | ----------: | :-----------: | :------------: |
|     0.1 |       459.6 |       406.4 |       no      |        no      |
|     0.5 |     1 092.9 |       828.4 |       no      |        no      |
|     1.0 |     1 873.3 |     1 348.4 |       no      |        no      |
|     2.0 |     3 398.2 |     2 364.1 |       **yes** |        no      |
|     5.0 |     7 721.5 |     5 237.1 |       **yes** |       **yes**  |

Interpretation:

* `Te` always >> `Tl` during the pulse — confirms the electrons absorb
  the laser energy first.
* In the prototype, the **lattice** crosses `T_m = 1687 K` somewhere
  between Ep = 1 µJ and 2 µJ, and crosses `T_v = 3538 K` between
  Ep = 2 µJ and 5 µJ.
* The single-temperature analytic estimate predicts melt at 0.42 µJ
  and vapor at 0.97 µJ — i.e. the TTM raises the apparent `Ep_melt`
  and `Ep_vap` by roughly **4×**.  This factor is realistic in
  principle (energy stored in the hot electron sub-system never reaches
  the lattice peak) but its exact value depends on the placeholder G,
  Ce and ke.

See:

* `results/v16_30ps/figures/v16_ttm_surface_temperature.png`
* `results/v16_30ps/figures/v16_ttm_depth_profiles.png`
* `results/v16_30ps/figures/v16_ttm_threshold_map.png`

---

## 6. Why the LS-PrePost view shows the hot zone in the upper-left

The V1 / V1.5 / V2 / V3A models are all **2-D axisymmetric**:

* `r = 0` is the symmetry axis and is drawn as the **left** edge.
* The wafer surface is drawn as the **top** edge (`z = THICKNESS`).
* The laser hits the wafer on-axis at the top surface, i.e. at the
  **corner of the rendered rectangle**.

In LS-PrePost's default view that corner shows up at the upper-left.
It is **not** "off-centre heating" — it's the symmetry axis of a
disk seen from the side.  Mirroring across `r = 0` puts the hot zone
in the **top centre** of a full cross-section, and revolving the
half-section about the z axis puts it in the **centre of the disk's
top face**.

`results/v16_30ps/figures/v16_axisymmetric_coordinate_explanation.png`
shows the three views side by side using the actual V1.5 / V2
`vapor_confirm_100ns` peak snapshot.  A zoomed coordinate-explanation
figure is also provided to make the hot zone near the axis/top surface
visually clearer:
`results/v16_30ps/figures/v16_axisymmetric_coordinate_explanation_zoomed.png`.

---

## 7. Files added by V1.6

| Layer       | Path |
| ----------- | --- |
| Config      | `config/v16_30ps.toml` |
| Script 1    | `scripts/estimate_v16_30ps_scales.py` |
| Script 2    | `scripts/ttm_1d_30ps_prototype.py` |
| Script 3    | `scripts/explain_axisymmetric_center_v16.py` |
| CSV 1       | `results/v16_30ps/v16_30ps_scale_summary.csv` |
| CSV 2       | `results/v16_30ps/v16_30ps_threshold_summary.csv` |
| CSV 3       | `results/v16_30ps/ttm_surface_temperature.csv` |
| CSV 4       | `results/v16_30ps/ttm_depth_profiles.csv` |
| Figure 1    | `results/v16_30ps/figures/v16_diffusion_length_vs_tau.png` |
| Figure 2    | `results/v16_30ps/figures/v16_threshold_energy_vs_tau.png` |
| Figure 3    | `results/v16_30ps/figures/v16_mesh_resolution_warning.png` |
| Figure 4    | `results/v16_30ps/figures/v16_ttm_surface_temperature.png` |
| Figure 5    | `results/v16_30ps/figures/v16_ttm_depth_profiles.png` |
| Figure 6    | `results/v16_30ps/figures/v16_ttm_threshold_map.png` |
| Figure 7    | `results/v16_30ps/figures/v16_axisymmetric_coordinate_explanation.png` |
| Figure 7b   | `results/v16_30ps/figures/v16_axisymmetric_coordinate_explanation_zoomed.png` |
| Doc         | `docs/README_v16_30ps.md` (this file) |

---

## 8. How to run

```powershell
.\.venv\Scripts\python.exe scripts\estimate_v16_30ps_scales.py
.\.venv\Scripts\python.exe scripts\ttm_1d_30ps_prototype.py
.\.venv\Scripts\python.exe scripts\explain_axisymmetric_center_v16.py
```

The first two are independent; the third reads V1.5 results read-only
(if `results/v15/vapor_confirm_100ns/tprint` is missing it falls back
to a synthetic Gaussian and the figure is clearly labelled as such).

---

## 9. Should we proceed to V1.7?

**Yes**, with the following constraints determined by V1.6:

| V1.7 parameter        | recommendation (from V1.6)                                   |
| --------------------- | ------------------------------------------------------------ |
| geometry              | local: r ∈ [0, 100 – 150 µm], z ∈ [0, 5 – 10 µm]              |
| surface dz            | **5 – 20 nm** (target L_diff / 5 ≈ 10 nm)                    |
| depth-grading         | grow to ~100 nm by z = 1 µm, ~500 nm by z = 5 – 10 µm        |
| radial-grading        | dr_min = 0.5 – 1 µm at r = 0, dr_max ≤ 5 µm beyond r = 50 µm |
| time-window           | 0 – 5 ns (or shorter; bulk of action by 1 ns)                |
| Ep cases (LS-DYNA)    | start with 0.5, 1.0, 2.0, 5.0 µJ; **TTM raises threshold ~4×**, so plan to also try 2 – 10 µJ |
| physics model         | single-temperature LS-DYNA thermal (V1.7 limitation; document) |

V1.7 will be a **single-temperature** lattice model (LS-DYNA has no
built-in TTM).  The V1.6 TTM prototype provides the cross-check:
V1.7 lattice temperatures should agree with the V1.6 prototype's `Tl`
to within the same factor that the prototype G / Ce / ke uncertainty
allows.

---

## 10. What V1.6 deliberately does NOT do

* Does not run LS-DYNA.
* Does not introduce element erosion, SPH, plume dynamics, or
  structured light.
* Does not change V1 / V1.5 / V2 / V3A inputs, results, or figures.
* Does not pretend that the TTM prototype is a quantitative model
  for silicon at any specific wavelength — all flags on figures and
  CSVs carry `PROTOTYPE PARAMETERS`.
* Does not commit to a V1.7 mesh or Ep list automatically — that is
  determined by §9 *under user approval* before V1.7 begins.
