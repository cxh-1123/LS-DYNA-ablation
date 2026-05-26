"""
build_v3b_plume_model.py
========================

V3B / V3B.1 -- build parametric plume / shock / ejecta model from V2.6 V3B driver.

Reads:
  - config/plume_model_v3b_30ps.toml
  - results/v26_30ps_threshold_ablation/v3b_driver/v3b_driver_ep5p0uj.json

Writes:
  - results/v3b_30ps_plume/v3b_plume_shock_metrics.csv
  - results/v3b_30ps_plume/v3b_driver_used.json
  - results/v3b_30ps_plume/v3b_t0.txt

Pure Python post-processor.  Does NOT modify V1.7 / V2.6 results.

Run:
    .\\.venv\\Scripts\\python.exe scripts\\build_v3b_plume_model.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tomllib
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass


def load_driver(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(
            f"V3B driver JSON not found: {path}\n"
            "Run scripts/export_v3b_driver_from_v26.py first."
        )
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_driver_geometry(driver: dict) -> tuple[float, float, float]:
    """Return (R0_um, d0_um, t0_ps) from driver + recommendations."""
    R0 = float(driver["recommended_crater_radius_um"])
    d0 = float(driver["recommended_crater_depth_um"])
    t0 = float(driver["recommended_t0_ps"])
    return R0, d0, t0


def evaluate_row(
    time_ps: float,
    t0_ps: float,
    R0: float,
    d0: float,
    v0: float,
    b: float,
    c0: float,
    beta: float,
    shock_ahead_min: float,
    density_decay_power: float,
    density_decay_ns: float,
    sigma_r_frac: float,
    sigma_z_frac: float,
    z_center_frac: float,
    min_sigma_r: float,
    min_sigma_z: float,
    ejecta_enabled: bool,
    max_ejecta_h: float,
    radial_spread: float,
    ejecta_decay_ns: float,
) -> dict:
    time_ns = time_ps * 1e-3
    dt_ns = max((time_ps - t0_ps) * 1e-3, 0.0)
    plume_active = dt_ns > 0.0

    if not plume_active:
        R_plume = R0
        R_shock = R0
        sigma_r = min_sigma_r
        sigma_z = min_sigma_z
        z_center = 0.0
        n0_rel = 0.0
        ejecta_h = 0.0
        ejecta_r = R0
        ejecta_active = False
        shock_ahead = 0.0
    else:
        R_plume = R0 + v0 * (dt_ns ** b)
        R_shock_raw = R0 + c0 * dt_ns + beta * np.sqrt(dt_ns)
        R_shock = max(R_shock_raw, R_plume + shock_ahead_min)
        shock_ahead = R_shock - R_plume

        sigma_r = max(min_sigma_r, sigma_r_frac * R_plume)
        sigma_z = max(min_sigma_z, sigma_z_frac * R_plume)
        z_center = z_center_frac * R_plume
        n0_rel = 1.0 / ((1.0 + dt_ns / density_decay_ns) ** density_decay_power)

        ejecta_active = ejecta_enabled
        if ejecta_enabled:
            ejecta_h = (
                max_ejecta_h
                * (1.0 - np.exp(-dt_ns / 0.2))
                * np.exp(-dt_ns / ejecta_decay_ns)
            )
            ejecta_r = R0 + radial_spread * (1.0 - np.exp(-dt_ns / 1.0))
        else:
            ejecta_h = 0.0
            ejecta_r = R0
            ejecta_active = False

    return {
        "time_ps": float(time_ps),
        "time_ns": float(time_ns),
        "dt_after_t0_ns": float(dt_ns),
        "plume_active": int(plume_active),
        "ejecta_active": int(ejecta_active),
        "R0_um": float(R0),
        "crater_depth_um": float(d0),
        "R_plume_um": float(R_plume),
        "R_shock_um": float(R_shock),
        "shock_ahead_um": float(shock_ahead),
        "sigma_r_um": float(sigma_r),
        "sigma_z_um": float(sigma_z),
        "z_center_um": float(z_center),
        "n0_rel": float(n0_rel),
        "ejecta_height_um": float(ejecta_h),
        "ejecta_radius_um": float(ejecta_r),
    }


def main() -> int:
    project_root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--toml", type=Path,
        default=project_root / "config" / "plume_model_v3b_30ps.toml",
    )
    ap.add_argument(
        "--out-root", type=Path,
        default=project_root / "results" / "v3b_30ps_plume",
    )
    args = ap.parse_args()

    if not args.toml.is_file():
        print(f"[ERROR] config not found: {args.toml}", file=sys.stderr)
        return 2

    with args.toml.open("rb") as fh:
        cfg = tomllib.load(fh)

    driver_path = project_root / cfg["driver"]["driver_json"]
    driver = load_driver(driver_path)
    R0, d0, t0_ps = resolve_driver_geometry(driver)

    t_end_ps = float(cfg["time"]["t_end_ns"]) * 1000.0
    dt_ps = float(cfg["time"]["dt_ps"])
    n_steps = int(round((t_end_ps - float(cfg["time"]["t_start_ps"])) / dt_ps)) + 1
    times_ps = np.linspace(float(cfg["time"]["t_start_ps"]), t_end_ps, n_steps)

    plume = cfg["plume"]
    shock = cfg["shock"]
    ejecta = cfg["ejecta"]
    density_decay_ns = float(plume.get("density_decay_ns", 2.0))

    rows = [
        evaluate_row(
            float(t), t0_ps, R0, d0,
            v0=float(plume["v0_um_per_ns_power"]),
            b=float(plume["b"]),
            c0=float(shock["c0_um_per_ns"]),
            beta=float(shock["beta_um_per_sqrt_ns"]),
            shock_ahead_min=float(shock["shock_ahead_min_um"]),
            density_decay_power=float(plume["density_decay_power"]),
            density_decay_ns=density_decay_ns,
            sigma_r_frac=float(plume["sigma_r_fraction"]),
            sigma_z_frac=float(plume["sigma_z_fraction"]),
            z_center_frac=float(plume["z_center_fraction"]),
            min_sigma_r=float(plume["min_sigma_r_um"]),
            min_sigma_z=float(plume["min_sigma_z_um"]),
            ejecta_enabled=bool(ejecta["enabled"]),
            max_ejecta_h=float(ejecta["max_ejecta_height_um"]),
            radial_spread=float(ejecta["radial_spread_um"]),
            ejecta_decay_ns=float(ejecta["decay_ns"]),
        )
        for t in times_ps
    ]

    args.out_root.mkdir(parents=True, exist_ok=True)
    out_csv = args.out_root / "v3b_plume_shock_metrics.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    used_driver = dict(driver)
    used_driver["_v3b_geometry_R0_um"] = R0
    used_driver["_v3b_geometry_d0_um"] = d0
    used_driver["_v3b_t0_ps"] = t0_ps
    used_driver["_v3b_version"] = cfg.get("meta", {}).get("version", "V3B")
    driver_used_path = args.out_root / "v3b_driver_used.json"
    driver_used_path.write_text(
        json.dumps(used_driver, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    t0_path = args.out_root / "v3b_t0.txt"
    t0_path.write_text(
        f"t0_ps={t0_ps:.6f}\nt0_ns={t0_ps * 1e-3:.9f}\n",
        encoding="utf-8",
    )

    R_plume = np.array([r["R_plume_um"] for r in rows])
    R_shock = np.array([r["R_shock_um"] for r in rows])
    ejecta_h = np.array([r["ejecta_height_um"] for r in rows])
    times_ns = np.array([r["time_ns"] for r in rows])

    version = cfg.get("meta", {}).get("version", "V3B")
    print("=" * 78)
    print(f"{version} -- build plume / shock / ejecta parametric model (30 ps driver)")
    print(f"  driver        : {driver_path}")
    print(f"  case          : {driver.get('case_name')}")
    print(f"  t0_ps         : {t0_ps:.3f}  (max crater-depth frame)")
    print(f"  R0_um         : {R0:.4f}")
    print(f"  d0_um         : {d0:.4f}")
    print(f"  output CSV    : {out_csv}")
    print(f"  time grid     : 0 -- {t_end_ps:.0f} ps, dt={dt_ps} ps ({n_steps} rows)")
    print("=" * 78)
    print(f"[OK] {out_csv}")
    print(f"[OK] {driver_used_path}")
    print(f"[OK] {t0_path}")

    print("\nR_plume / R_shock snapshots:")
    for label, t_ns in [
        ("0.1 ns", 0.1), ("0.5 ns", 0.5), ("1 ns", 1.0),
        ("2 ns", 2.0), ("5 ns", 5.0),
    ]:
        i = int(np.argmin(np.abs(times_ns - t_ns)))
        ahead = R_shock[i] - R_plume[i]
        print(
            f"  t={label:>6}  R_plume={R_plume[i]:7.2f} um  "
            f"R_shock={R_shock[i]:7.2f} um  ahead={ahead:.2f} um"
        )
    print(f"\n  max ejecta_height_um : {float(ejecta_h.max()):.2f} um")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
