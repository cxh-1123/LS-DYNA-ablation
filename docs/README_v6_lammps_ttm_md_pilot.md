# V6 — LAMMPS TTM-MD pilot (Si, 30 ps, ep_5p0uj)

V6 is the **first atomistic TTM+MD step** toward paper-like plume/ejecta and
structured-light (V4) inputs.  It does **not** modify V1.7 / V2.6 / V3B / V5A.

---

## 1. Position in the pipeline

```
V1.6 TTM-only (1D)     →  scale / Te-Tl timing
V1.7 LS-DYNA (single-T) →  bulk / ns reference
V2.6 / V3B             →  threshold + parametric proxy (until V6 data exists)
V5A                    →  LS-DYNA marker plume (ns, d3plot)
V6 (this)              →  LAMMPS TTM-MD ps plume / ejecta / density → V4
```

---

## 2. Aligned parameters (ep_5p0uj)

| Parameter | Value | Source |
| --- | ---: | --- |
| Ep | 5.0 µJ | V17 / V2.6 main case |
| τ pulse | 30 ps | V16 / V17 |
| w₀ | 35 µm | V16 / V17 |
| A | 0.5 | V17 |
| δ_abs (266 nm) | 5 nm | V16 placeholder |
| T_melt / T_vap | 1687 / 3538 K | project thresholds |

TTM electron parameters (`C_e`, `G`, `k_e`) follow **V1.6 placeholders** —
calibrate before quantitative claims.

---

## 3. Files

| Path | Role |
| --- | --- |
| `config/v6_lammps_ttm_md_pilot.toml` | Pilot knobs |
| `scripts/build_v6_lammps_ttm_input.py` | Generate LAMMPS input |
| `scripts/run_v6_selected.ps1` | Run LAMMPS |
| `scripts/plot_v6_lammps_ttm_snapshots.py` | Snapshot figures |
| `scripts/export_v6_density_for_v4.py` | V4 handoff CSV |
| `models/v6_lammps_ttm_md_pilot/ep_5p0uj/` | Generated inputs |
| `results/v6_lammps_ttm_md_pilot/ep_5p0uj/` | LAMMPS outputs |
| `potentials/Si.sw` | Stillinger–Weber (user-supplied) |

---

## 4. Build

```powershell
python scripts\build_v6_lammps_ttm_input.py
```

Generates:

- `in_ep_5p0uj.lammps`
- `laser_ttm_source_ep_5p0uj.txt`
- `v6_run_manifest.json`

---

## 5. Prerequisites

1. **LAMMPS** with `fix ttm` (MOLECULE package)
2. **`potentials/Si.sw`** — see `potentials/README.md`
3. Enough CPU/RAM for ~10⁴–10⁵ atoms (pilot slab)

Check:

```powershell
lmp -help | findstr ttm
```

---

## 6. Run

Dry-run:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v6_selected.ps1 -DryRun
```

Execute:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_v6_selected.ps1 -OnlyCase ep_5p0uj
```

Or manually from `models/v6_lammps_ttm_md_pilot/ep_5p0uj/`:

```powershell
lmp -in in_ep_5p0uj.lammps
```

Snapshots: `snapshot_0ps.lammpstrj` … `snapshot_500ps.lammpstrj`

---

## 7. Post-process

```powershell
python scripts\plot_v6_lammps_ttm_snapshots.py --case ep_5p0uj
python scripts\export_v6_density_for_v4.py --case ep_5p0uj
```

V4 will read:

- `results/v6_lammps_ttm_md_pilot/ep_5p0uj/v6_density_grid_*ps.csv`
- `v6_v4_handoff_manifest.json`

Until V6 finishes, V4 can still use **V3B proxy** (`v3b_plume_shock_metrics.csv`).

---

## 8. Unified figure handoff (V5B frame)

Target snapshot times match V5B main sequence:

`0, 10, 30, 50, 100, 200, 500 ps`

Fixed display frame (for comparison plots):

- r: [-60, 60] µm  
- z: [-8, 80] µm  

After V6 snapshots exist, add a converter script (V6.1) to plot on the **same**
layout as `v5b_early_ablation_plume_sequence.png`.

---

## 9. Limitations (pilot)

- Small MD box (~18×18×36 diamond cells) — not full 500 µm wafer
- TTM parameters are **qualitative** (V1.6 order)
- Spatial laser source at **peak**; full 30 ps Gaussian time envelope → V6.1
- No automatic coupling back to LS-DYNA `.k` files
- **Not** a replacement for peer-review TTM-MD until Si.sw + G + Ce calibrated

---

## 10. Roadmap

| Step | Content |
| --- | --- |
| **V6.0** (now) | Build + run pilot, export density for V4 |
| **V6.1** | Time-dependent TTM source, larger box, melt/ vapor analysis |
| **V6.2** | Unified V5B-layout plot from real MD snapshots |
| **V4** | Structured light using V6 density (or V3B fallback) |

---

## 11. What stays unchanged

- V1.7 / V2.6 extract logic and CSV values  
- V3B / V5A scripts and results (unless you explicitly re-run them)  
- LS-DYNA `.k` / tprint archives  
