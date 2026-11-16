"""Refit the adopted OAK progression mixtures with the death margins locked.

This is a deterministic curve fit, not a Monte Carlo trial simulation.  The
formal reference parameters are used only after optimization for comparison.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
from dataclasses import asdict
from pathlib import Path
import sys
_SUPPORT_CODE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SUPPORT_CODE_ROOT / 'src/common'))
from output_paths import require_external_output
import sys
_SOURCE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SOURCE_ROOT / 'src/dgm'))

import numpy as np
import scipy
from scipy.optimize import least_squares

import oak_fit_core as core


def unpack(theta: np.ndarray, rho: float, death: dict) -> core.MixtureParameters:
    pi, mu1, log_sigma1, log_delta, log_sigma2 = map(float, theta)
    return core.MixtureParameters(
        rho=rho,
        pi=pi,
        mu_r1=mu1,
        sigma_r1=math.exp(log_sigma1),
        mu_r2=mu1 + math.exp(log_delta),
        sigma_r2=math.exp(log_sigma2),
        mu_d=float(death["mu_d"]),
        sigma_d=float(death["sigma_d"]),
    )


def optimizer_record(result, initial: list, candidate: str) -> dict:
    return {
        "candidate": candidate,
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nfev": int(result.nfev),
        "sum_squared_residuals": float(np.dot(result.fun, result.fun)),
        "initial_theta": list(map(float, initial)),
        "optimized_theta": list(map(float, result.x)),
    }


def fit_control(curve, death: dict, settings: dict, rho: float):
    risk_t = np.asarray(settings["risk_times"], dtype=float)
    risk_n = np.asarray(settings["risk_counts"], dtype=float)
    grid = np.linspace(0.001, 18.0, 360)
    target = core.interp_curve(curve, grid)
    weights = np.sqrt(np.maximum(np.interp(grid, risk_t, risk_n / 425.0), 0.1))
    checkpoints = np.array([6.0, 12.0, 18.0])
    checkpoint_target = core.interp_curve(curve, checkpoints)
    median_target = float(settings["median_anchor"])
    sigma_floor = float(settings["sigma_lower_bound"])

    def residual(theta):
        p = unpack(theta, rho, death)
        survival = core.mixture_pfs_survival(grid, p)
        anchors = core.mixture_pfs_survival(np.r_[median_target, checkpoints], p)
        return np.r_[
            (survival - target) * weights / 0.012,
            np.repeat((anchors[0] - 0.5) / 0.008, 50),
            np.repeat((anchors[1:] - checkpoint_target) / 0.015, 12),
        ]

    lower = [0.05, math.log(0.1), math.log(sigma_floor), math.log(0.05), math.log(sigma_floor)]
    upper = [0.95, math.log(20), math.log(3), math.log(30), math.log(3)]
    initial = [0.11, 0.33, math.log(max(0.1, sigma_floor)), math.log(1.59 - 0.33), math.log(0.88)]
    result = least_squares(
        residual, initial, bounds=(lower, upper), max_nfev=3000,
        xtol=1e-11, ftol=1e-11, gtol=1e-11,
    )
    if not result.success:
        raise RuntimeError(f"Control optimization failed: {result.message}")
    return unpack(result.x, rho, death), optimizer_record(result, initial, "control_sigma_floor_0.30")


def fit_treatment(curve, death: dict, settings: dict, rho: float, control):
    risk_t = np.asarray(settings["risk_times"], dtype=float)
    risk_n = np.asarray(settings["risk_counts"], dtype=float)
    grid = np.linspace(0.001, 18.0, 320)
    target = core.interp_curve(curve, grid)
    weights = np.sqrt(np.maximum(np.interp(grid, risk_t, risk_n / risk_n[0]), 0.1))
    checkpoints = np.array([6.0, 12.0, 18.0])
    checkpoint_target = core.interp_curve(curve, checkpoints)
    hr_grid = np.linspace(0.0, 24.0, 1601)
    control_survival = core.mixture_pfs_survival(hr_grid, control)
    median_target = float(settings["median_anchor"])
    hr_target = float(settings["hr_optimization_target"])

    def residual(theta):
        p = unpack(theta, rho, death)
        survival = core.mixture_pfs_survival(grid, p)
        anchors = core.mixture_pfs_survival(np.r_[median_target, checkpoints], p)
        hr = core.cox_hr_from_survival(
            hr_grid, core.mixture_pfs_survival(hr_grid, p), control_survival,
        )
        return np.r_[
            (survival - target) * weights / 0.018,
            np.repeat((anchors[0] - 0.5) / 0.015, 30),
            np.repeat((anchors[1:] - checkpoint_target) / 0.025, 8),
            np.repeat(math.log(hr / hr_target) / 0.02, 45),
        ]

    lower = [0.05, math.log(0.08), math.log(0.05), math.log(0.05), math.log(0.05)]
    upper = [0.95, math.log(12), math.log(3), math.log(30), math.log(3)]
    # The documented first start was the earlier fit, not the final reference.
    initial = [
        0.2264164067590971,
        0.3535388761562825,
        math.log(0.08489730075749218),
        math.log(1.5295924778235532 - 0.3535388761562825),
        math.log(1.3068570883563746),
    ]
    result = least_squares(
        residual, initial, bounds=(lower, upper), max_nfev=3500,
        xtol=1e-11, ftol=1e-11, gtol=1e-11,
    )
    if not result.success:
        raise RuntimeError(f"Treatment optimization failed: {result.message}")
    return unpack(result.x, rho, death), optimizer_record(result, initial, "treatment_first_local_candidate")


def read_reference(path: Path) -> dict[str, core.MixtureParameters]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    result = {}
    for arm in ["control", "treatment"]:
        components = [r for r in rows if r["case"] == "OAK" and r["scenario"] == "H1-observed" and r["arm"] == arm]
        components.sort(key=lambda x: int(x["component"]))
        if len(components) != 2:
            raise ValueError(f"Expected two OAK H1-observed reference components for {arm}")
        c1, c2 = components
        for field in ["mu_log_D", "sigma_log_D", "rho_log_R_D"]:
            if float(c1[field]) != float(c2[field]):
                raise ValueError(f"The OAK reference does not share {field} across components")
        if abs(float(c1["weight"]) + float(c2["weight"]) - 1.0) > 1e-12:
            raise ValueError("Reference mixture weights do not sum to 1")
        result[arm] = core.MixtureParameters(
            float(c1["rho_log_R_D"]), float(c1["weight"]),
            float(c1["mu_log_R"]), float(c1["sigma_log_R"]),
            float(c2["mu_log_R"]), float(c2["sigma_log_R"]),
            float(c1["mu_log_D"]), float(c1["sigma_log_D"]),
        )
    return result


def diagnostics(curve, fitted, control, settings: dict) -> dict:
    grid = np.linspace(0.001, 18.0, 6000)
    survival = core.mixture_pfs_survival(grid, fitted)
    median = core.survival_median(grid, survival)
    checkpoints = np.array([6.0, 12.0, 18.0])
    errors = np.abs(core.mixture_pfs_survival(checkpoints, fitted) - core.interp_curve(curve, checkpoints))
    result = {
        "pfs_median_months": median,
        "median_relative_error": abs(median / float(settings["median_anchor"]) - 1.0),
        "pfs_max_absolute_curve_error_through_18_months": float(np.max(np.abs(survival - core.interp_curve(curve, grid)))),
        "pfs_absolute_error_at_6_12_18_months": list(map(float, errors)),
    }
    if control is not None:
        hr_grid = np.linspace(0.0, 24.0, 5001)
        hr = core.cox_hr_from_survival(hr_grid, core.mixture_pfs_survival(hr_grid, fitted), core.mixture_pfs_survival(hr_grid, control))
        result["curve_level_cox_hr_through_24_months"] = hr
        result["hr_relative_error_to_ipd_anchor"] = abs(hr / float(settings["hr_validation_anchor"]) - 1.0)
    result["historical_curve_acceptance_pass"] = bool(
        result["median_relative_error"] <= 0.05
        and result["pfs_max_absolute_curve_error_through_18_months"] <= 0.05
        and max(result["pfs_absolute_error_at_6_12_18_months"]) <= 0.05
        and result.get("hr_relative_error_to_ipd_anchor", 0.0) <= 0.05
    )
    return result


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir()) / "config/dgm_fitting/oak")
    parser.add_argument("--curve-dir", type=Path, default=next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir()) / "inputs/clinical_data/OAK", help="Directory containing the existing OAK digitized curves.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for newly computed results; must be outside the package and input directories.")
    parser.add_argument("--parameter-atol", type=float, default=5e-5, help="Absolute parameter tolerance for numerical refitting versus the archived reference.")
    parser.add_argument("--curve-atol", type=float, default=1e-5, help="Maximum absolute fitted-reference survival difference.")
    args = parser.parse_args()
    package_root = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
    inputs = args.input_dir.resolve()
    curve_dir = args.curve_dir.resolve()
    output = require_external_output(args.output_dir, _SUPPORT_CODE_ROOT)
    for protected in [package_root, inputs, curve_dir]:
        if output == protected or output.is_relative_to(protected):
            parser.error("--output-dir must be outside the package and input directories")
    if not all(math.isfinite(tolerance) and tolerance > 0 for tolerance in [args.parameter_atol, args.curve_atol]):
        parser.error("Comparison tolerances must be finite and positive")
    settings = json.loads((inputs / "fit_settings.json").read_text(encoding="utf-8"))
    death_lock = json.loads((inputs / "oak_locked_death_parameters.json").read_text(encoding="utf-8"))
    curves = {arm: core.read_curve(curve_dir / f"OAK_A_ITT_pfs_{arm}.csv") for arm in ["control", "treatment"]}
    rho = float(settings["rho"])
    control, control_optimizer = fit_control(curves["control"], death_lock["control"], settings["control"], rho)
    treatment, treatment_optimizer = fit_treatment(curves["treatment"], death_lock["treatment"], settings["treatment"], rho, control)
    fitted = {"control": control, "treatment": treatment}
    reference = read_reference(inputs / "locked_reference.csv")
    comparison = []
    curve_comparison = {}
    grid = np.linspace(0.0, 27.0, 1001)
    curve_rows = []
    for arm, p in fitted.items():
        for name, value in asdict(p).items():
            expected = getattr(reference[arm], name)
            difference = abs(value - expected)
            comparison.append({"arm": arm, "parameter": name, "fitted": value, "locked_reference": expected, "absolute_difference": difference, "within_tolerance": difference <= args.parameter_atol})
        pfs = core.mixture_pfs_survival(grid, p)
        pfs_reference = core.mixture_pfs_survival(grid, reference[arm])
        os = core.lognormal_survival(grid, p.mu_d, p.sigma_d)
        os_reference = core.lognormal_survival(grid, reference[arm].mu_d, reference[arm].sigma_d)
        curve_comparison[arm] = {"pfs_max_absolute_difference": float(np.max(np.abs(pfs-pfs_reference))), "os_max_absolute_difference": float(np.max(np.abs(os-os_reference)))}
        for i, time in enumerate(grid):
            curve_rows.append({"arm": arm, "time_months": time, "fitted_pfs": pfs[i], "reference_pfs": pfs_reference[i], "fitted_os": os[i], "reference_os": os_reference[i]})
    checks = {"control": diagnostics(curves["control"], control, None, settings["control"]), "treatment": diagnostics(curves["treatment"], treatment, control, settings["treatment"])}
    passed = all(r["within_tolerance"] for r in comparison) and all(v <= args.curve_atol for arm in curve_comparison.values() for v in arm.values()) and all(r["historical_curve_acceptance_pass"] for r in checks.values())
    report = {
        "case": "OAK", "scenario": "H1-observed", "status": "PASS" if passed else "FAIL",
        "method": "OS margins fixed; arm-specific progression-mixture parameters fitted to PFS",
        "free_parameters_per_arm": ["pi", "mu_R1", "log_sigma_R1", "log_delta_mu_R", "log_sigma_R2"],
        "death_parameters_optimized": False, "os_residual_in_objective": False,
        "candidate_adoption": settings["candidate_adoption"],
        "optimizer": {"control": control_optimizer, "treatment": treatment_optimizer},
        "fitted_parameters": {arm: asdict(p) for arm, p in fitted.items()},
        "curve_diagnostics": checks, "fitted_reference_curve_differences": curve_comparison,
        "comparison_tolerances": {"parameter_absolute": args.parameter_atol, "curve_absolute": args.curve_atol},
        "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "scope": "Deterministic refit and reference check only; no study simulation, MCMC, or locked configuration is changed.",
    }
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "fitted_parameters.csv", [{"arm": arm, **asdict(p)} for arm,p in fitted.items()])
    write_csv(output / "parameter_comparison.csv", comparison)
    write_csv(output / "fitted_reference_curves.csv", curve_rows)
    (output / "fit_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"case": "OAK", "status": report["status"], "control_pi": control.pi, "treatment_pi": treatment.pi, "report": str(output / "fit_report.json")}, indent=2))
    return 0 if passed else 2


if __name__ == "__main__":
    sys.exit(main())
