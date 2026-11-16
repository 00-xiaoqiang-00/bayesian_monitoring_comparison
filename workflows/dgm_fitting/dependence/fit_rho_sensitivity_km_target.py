"""Formal rho-sensitivity refit against original/digitized observable curves.

The fitting target is never a curve generated from the locked DGM.  For each
case and treatment arm, the archived H1-observed curve is used directly:

* PFS for IMpower131, OAK, and KEYNOTE-598;
* recurrence cumulative incidence for survival_colon;
* OS for all four cases.

The observed control/treatment DGM is refit at baseline, low, and high rho
with the archived death marginal and mixture weights fixed.  H0 and
H1-standard are then derived from the refitted observed control arm using the
locked standard-effect contrasts from E07.  Neither scenario is refit to an
artificial DGM-generated curve, and no decision threshold is recalibrated here.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import sys
_SUPPORT_CODE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SUPPORT_CODE_ROOT / 'src/common'))
from output_paths import require_external_output

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy.special import ndtr


REPO_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
LOCK = REPO_ROOT / "config/dgm_fitting/dependence/settings.json"
DGM_INPUT = REPO_ROOT / "config/dgm_parameters.csv"
TARGET_DIR = None
TARGET_MANIFEST = None
OUT = None
SELECTED_CASES = None


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def lognormal_surv(t: float | np.ndarray, mu: float, sigma: float) -> np.ndarray:
    t_arr = np.asarray(t, dtype=float)
    out = np.ones_like(t_arr)
    positive = t_arr > 0
    z = (np.log(np.maximum(t_arr[positive], 1e-12)) - mu) / sigma
    out[positive] = ndtr(-z)
    return out


def bivar_pfs_surv(t: float, mu_r: float, sig_r: float, mu_d: float, sig_d: float, rho: float) -> float:
    """P(R>t,D>t) using one-dimensional conditional integration."""
    if t <= 0:
        return 1.0
    z_d = (math.log(t) - mu_d) / sig_d
    lower = max(z_d, -12.0)

    def integrand(z: float) -> float:
        cond_mean = mu_r + rho * sig_r * z
        cond_sd = sig_r * math.sqrt(max(1e-12, 1.0 - rho * rho))
        p_r = float(ndtr(-(math.log(t) - cond_mean) / cond_sd))
        return math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi) * p_r

    return float(quad(integrand, lower, 12.0, epsabs=3e-8, epsrel=3e-7, limit=120)[0])


def recurrence_cif(t: float, mu_r: float, sig_r: float, mu_d: float, sig_d: float, rho: float) -> float:
    """P(R<=t, R<D) for a bivariate lognormal semi-competing-risk model."""
    if t <= 0:
        return 0.0
    upper = min(12.0, (math.log(t) - mu_r) / sig_r)
    if upper <= -12.0:
        return 0.0

    def integrand(z: float) -> float:
        log_r = mu_r + sig_r * z
        cond_mean = mu_d + rho * sig_d * z
        cond_sd = sig_d * math.sqrt(max(1e-12, 1.0 - rho * rho))
        p_d_after = float(ndtr(-(log_r - cond_mean) / cond_sd))
        return math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi) * p_d_after

    return float(quad(integrand, -12.0, upper, epsabs=3e-8, epsrel=3e-7, limit=120)[0])


def model_curve(times: np.ndarray, components: list[dict[str, float]], rho: float, endpoint: str) -> np.ndarray:
    out = np.zeros(len(times), dtype=float)
    if endpoint == "OS":
        for comp in components:
            out += comp["weight"] * lognormal_surv(times, comp["mu_D"], comp["sigma_D"])
        return np.clip(out, 0.0, 1.0)
    if endpoint == "recurrence_cif":
        if len(components) != 1:
            raise ValueError("recurrence_cif target requires one shared latent component")
        comp = components[0]
        return np.asarray([
            recurrence_cif(t, comp["mu_R"], comp["sigma_R"], comp["mu_D"], comp["sigma_D"], rho)
            for t in times
        ])
    for comp in components:
        out += comp["weight"] * np.asarray([
            bivar_pfs_surv(t, comp["mu_R"], comp["sigma_R"], comp["mu_D"], comp["sigma_D"], rho)
            for t in times
        ])
    return np.clip(out, 0.0, 1.0)


def target_points(rows: list[dict[str, str]], max_points: int = 220) -> tuple[np.ndarray, np.ndarray]:
    times = np.asarray([float(r["time_months"]) for r in rows], dtype=float)
    values = np.asarray([float(r["curve_value"]) for r in rows], dtype=float)
    if len(times) > max_points:
        idx = np.unique(np.round(np.linspace(0, len(times) - 1, max_points)).astype(int))
        return times[idx], values[idx]
    return times, values


def component_rows(rows: list[dict[str, str]]) -> list[dict[str, float]]:
    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[int(row["component"])].append(row)
    output = []
    for component in sorted(grouped):
        row = grouped[component][0]
        output.append({
            "component": float(component),
            "weight": as_float(row, "weight"),
            "mu_R": as_float(row, "mu_log_R"),
            "sigma_R": as_float(row, "sigma_log_R"),
            "mu_D": as_float(row, "mu_log_D"),
            "sigma_D": as_float(row, "sigma_log_D"),
        })
    return output


def fit_progression(target_t: np.ndarray, target_y: np.ndarray, start: list[dict[str, float]], rho: float, endpoint: str) -> tuple[list[dict[str, float]], dict[str, object]]:
    x0 = []
    for comp in start:
        x0.extend([comp["mu_R"], math.log(comp["sigma_R"])])
    x0 = np.asarray(x0, dtype=float)

    def unpack(x: np.ndarray) -> list[dict[str, float]]:
        return [
            {
                "component": start[k]["component"],
                "weight": start[k]["weight"],
                "mu_R": float(x[2 * k]),
                "sigma_R": float(math.exp(x[2 * k + 1])),
                "mu_D": start[k]["mu_D"],
                "sigma_D": start[k]["sigma_D"],
            }
            for k in range(len(start))
        ]

    def residual(x: np.ndarray) -> np.ndarray:
        candidate = unpack(x)
        fitted = model_curve(target_t, candidate, rho, endpoint)
        curve = (fitted - target_y) / 0.01
        # Small stabilizer prevents a non-identifiable mixture from moving far
        # from the archived curve-fit representation when the curves are flat.
        regularization = []
        for k, comp in enumerate(start):
            regularization.extend([
                0.01 * (x[2 * k] - comp["mu_R"]),
                0.01 * (x[2 * k + 1] - math.log(comp["sigma_R"])),
            ])
        return np.concatenate([curve, np.asarray(regularization)])

    lower, upper = [], []
    for comp in start:
        lower.extend([comp["mu_R"] - 3.5, math.log(0.05)])
        upper.extend([comp["mu_R"] + 3.5, math.log(5.0)])
    result = least_squares(
        residual,
        x0=x0,
        bounds=(np.asarray(lower), np.asarray(upper)),
        max_nfev=2500,
        xtol=1e-10,
        ftol=1e-10,
        gtol=1e-10,
    )
    candidate = unpack(result.x)
    fitted = model_curve(target_t, candidate, rho, endpoint)
    delta = fitted - target_y
    return candidate, {
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_nfev": int(result.nfev),
        "optimizer_cost": float(result.cost),
        "weighted_rmse": float(math.sqrt(np.mean(delta * delta))),
        "mae": float(np.mean(np.abs(delta))),
        "max_abs": float(np.max(np.abs(delta))),
    }


def fit_death(target_t: np.ndarray, target_y: np.ndarray, start: list[dict[str, float]]) -> tuple[list[dict[str, float]], dict[str, object]]:
    """Fit the marginal death mixture once to the observed OS target.

    This is deliberately outside the rho loop.  The resulting death marginal
    is then held fixed for low/baseline/high rho, which separates observable
    OS fitting from the dependence sensitivity.
    """
    x0 = []
    for comp in start:
        x0.extend([comp["mu_D"], math.log(comp["sigma_D"])])
    x0 = np.asarray(x0, dtype=float)

    def unpack(x: np.ndarray) -> list[dict[str, float]]:
        return [
            {
                "component": start[k]["component"],
                "weight": start[k]["weight"],
                "mu_R": start[k]["mu_R"],
                "sigma_R": start[k]["sigma_R"],
                "mu_D": float(x[2 * k]),
                "sigma_D": float(math.exp(x[2 * k + 1])),
            }
            for k in range(len(start))
        ]

    def residual(x: np.ndarray) -> np.ndarray:
        candidate = unpack(x)
        fitted = model_curve(target_t, candidate, 0.0, "OS")
        curve = (fitted - target_y) / 0.01
        regularization = []
        for k, comp in enumerate(start):
            regularization.extend([
                0.01 * (x[2 * k] - comp["mu_D"]),
                0.01 * (x[2 * k + 1] - math.log(comp["sigma_D"])),
            ])
        return np.concatenate([curve, np.asarray(regularization)])

    lower, upper = [], []
    for comp in start:
        lower.extend([comp["mu_D"] - 3.5, math.log(0.05)])
        upper.extend([comp["mu_D"] + 3.5, math.log(5.0)])
    result = least_squares(
        residual,
        x0=x0,
        bounds=(np.asarray(lower), np.asarray(upper)),
        max_nfev=2500,
        xtol=1e-10,
        ftol=1e-10,
        gtol=1e-10,
    )
    candidate = unpack(result.x)
    fitted = model_curve(target_t, candidate, 0.0, "OS")
    delta = fitted - target_y
    return candidate, {
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_nfev": int(result.nfev),
        "optimizer_cost": float(result.cost),
        "weighted_rmse": float(math.sqrt(np.mean(delta * delta))),
        "mae": float(np.mean(np.abs(delta))),
        "max_abs": float(np.max(np.abs(delta))),
    }


def curve_metrics(target_t: np.ndarray, target_y: np.ndarray, fitted: np.ndarray) -> dict[str, float]:
    delta = fitted - target_y
    return {
        "weighted_rmse": float(math.sqrt(np.mean(delta * delta))),
        "mae": float(np.mean(np.abs(delta))),
        "max_abs": float(np.max(np.abs(delta))),
    }


def median_time(times: np.ndarray, values: np.ndarray) -> float:
    hit = np.flatnonzero(values <= 0.5)
    return float(times[hit[0]]) if hit.size else float("nan")


def row_to_component(row: dict[str, object]) -> dict[str, float]:
    return {
        "component": float(row["component"]),
        "weight": float(row["weight"]),
        "mu_R": float(row["mu_log_R"]),
        "sigma_R": float(row["sigma_log_R"]),
        "mu_D": float(row["mu_log_D"]),
        "sigma_D": float(row["sigma_log_D"]),
    }


def output_row(case: str, scenario: str, arm: str, rho_label: str, rho: float, comp: dict[str, float], derivation: str, source_scenario: str) -> dict[str, object]:
    return {
        "case": case,
        "scenario": scenario,
        "arm": arm,
        "rho_label": rho_label,
        "rho": rho,
        "component": int(comp["component"]),
        "weight": comp["weight"],
        "mu_log_R": comp["mu_R"],
        "median_R_months": math.exp(comp["mu_R"]),
        "sigma_log_R": comp["sigma_R"],
        "mu_log_D": comp["mu_D"],
        "median_D_months": math.exp(comp["mu_D"]),
        "sigma_log_D": comp["sigma_D"],
        "cov_log_R_D": rho * comp["sigma_R"] * comp["sigma_D"],
        "derivation": derivation,
        "source_scenario": source_scenario,
    }


def main() -> None:
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    if SELECTED_CASES:
        lock["rho_grid"] = {k:v for k,v in lock["rho_grid"].items() if k in SELECTED_CASES}
    target_manifest = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    target_rows = list(csv.DictReader((TARGET_DIR / "target_manifest.csv").open(encoding="utf-8"))) if (TARGET_DIR / "target_manifest.csv").exists() else target_manifest["target_rows"]
    # target_rows in target_manifest contains one row per curve and its file.
    targets: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for spec in target_rows:
        key = (spec["case"], spec["endpoint"], spec["arm"])
        path = Path(spec["target_file"])
        if not path.is_absolute():
            path = TARGET_DIR / path
        targets[key] = load_csv(path)

    dgm_rows = load_csv(DGM_INPUT)
    e09_rows = [r for r in dgm_rows if r["scenario"] == "H1-observed"]
    e07_rows = [r for r in dgm_rows if r["scenario"] == "H1-standard"]
    observed: dict[tuple[str, str], list[dict[str, float]]] = {}
    standard: dict[tuple[str, str], list[dict[str, float]]] = {}
    for case in lock["rho_grid"]:
        for arm in ("control", "treatment"):
            observed[(case, arm)] = component_rows([r for r in e09_rows if r["case"] == case and r["arm"] == arm])
            standard[(case, arm)] = component_rows([r for r in e07_rows if r["case"] == case and r["arm"] == arm])

    direct_parameter_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    direct_fits: dict[tuple[str, str, str], list[dict[str, float]]] = {}
    fixed_death: dict[tuple[str, str], list[dict[str, float]]] = {}
    death_check_metrics: dict[tuple[str, str], dict[str, object]] = {}
    # Death parameters are part of the locked observable curve-fit DGM.  They
    # are not re-estimated under rho sensitivity.  We only calculate their OS
    # discrepancy against the digitized/patient-level target as an audit check.
    for case in lock["rho_grid"]:
        for arm in ("control", "treatment"):
            target_os = targets[(case, "OS", arm)]
            target_os_t, target_os_y = target_points(target_os)
            fixed_death[(case, arm)] = observed[(case, arm)]
            fitted_os = model_curve(target_os_t, observed[(case, arm)], 0.0, "OS")
            metrics = curve_metrics(target_os_t, target_os_y, fitted_os)
            death_check_metrics[(case, arm)] = {
                **metrics,
                "optimizer_success": True,
                "optimizer_message": "locked E09 death marginal; no rho-specific death refit",
            }
    for case, grid in lock["rho_grid"].items():
        endpoint = lock["target_endpoints"][case][0]
        for arm in ("control", "treatment"):
            target_primary = targets[(case, endpoint, arm)]
            target_primary_t, target_primary_y = target_points(target_primary)
            target_os = targets[(case, "OS", arm)]
            target_os_t, target_os_y = target_points(target_os)
            starting_progression = [dict(comp) for comp in fixed_death[(case, arm)]]
            for rho_label in ("low", "baseline", "high"):
                rho = float(grid[rho_label])
                fitted, primary_metrics = fit_progression(
                    target_primary_t, target_primary_y, starting_progression, rho, endpoint
                )
                fitted_os = model_curve(target_os_t, fitted, rho, "OS")
                os_metrics = curve_metrics(target_os_t, target_os_y, fitted_os)
                direct_fits[(case, arm, rho_label)] = fitted
                for comp in fitted:
                    direct_parameter_rows.append(output_row(
                        case, "H1-observed", arm, rho_label, rho, comp,
                        "KM_TARGET_DIRECT_REFIT", "H1-observed",
                    ))
                metric_rows.append({
                    "case": case,
                    "scenario": "H1-observed",
                    "arm": arm,
                    "rho_label": rho_label,
                    "rho": rho,
                    "primary_endpoint": endpoint,
                    "primary_n_target_points": len(target_primary),
                    "primary_target_max_months": float(target_primary[-1]["time_months"]),
                    "primary_target_median_or_cif": median_time(target_primary_t, target_primary_y) if endpoint == "PFS" else float(target_primary_y[-1]),
                    "primary_fitted_median_or_cif": median_time(target_primary_t, model_curve(target_primary_t, fitted, rho, endpoint)) if endpoint == "PFS" else float(model_curve(target_primary_t, fitted, rho, endpoint)[-1]),
                    "primary_weighted_rmse": primary_metrics["weighted_rmse"],
                    "primary_mae": primary_metrics["mae"],
                    "primary_max_abs": primary_metrics["max_abs"],
                    "os_n_target_points": len(target_os),
                    "os_target_max_months": float(target_os[-1]["time_months"]),
                    "os_weighted_rmse": os_metrics["weighted_rmse"],
                    "os_mae": os_metrics["mae"],
                    "os_max_abs": os_metrics["max_abs"],
                    "death_check_weighted_rmse": death_check_metrics[(case, arm)]["weighted_rmse"],
                    "death_check_mae": death_check_metrics[(case, arm)]["mae"],
                    "death_check_max_abs": death_check_metrics[(case, arm)]["max_abs"],
                    "death_check_optimizer_success": death_check_metrics[(case, arm)]["optimizer_success"],
                    "optimizer_success": primary_metrics["optimizer_success"],
                    "optimizer_nfev": primary_metrics["optimizer_nfev"],
                    "optimizer_cost": primary_metrics["optimizer_cost"],
                    "optimizer_message": primary_metrics["optimizer_message"],
                    "target_source_kind": target_primary[0]["source_kind"],
                })

    # Derive H0 and H1-standard from the KM-fitted observed control arm. This
    # preserves the locked standard effect while avoiding a second artificial
    # curve-fitting target for H0/H1-standard.
    all_parameters = list(direct_parameter_rows)
    for case, grid in lock["rho_grid"].items():
        for rho_label in ("low", "baseline", "high"):
            rho = float(grid[rho_label])
            control = direct_fits[(case, "control", rho_label)]
            std_control = standard[(case, "control")]
            std_treatment = standard[(case, "treatment")]
            for arm in ("control", "treatment"):
                for comp in control:
                    all_parameters.append(output_row(
                        case, "H0", arm, rho_label, rho, dict(comp),
                        "DERIVED_FROM_KM_FITTED_CONTROL" if arm == "control" else "DERIVED_H0_TREATMENT_EQUALS_CONTROL",
                        "KM-fitted H1-observed control",
                    ))
            standard_treatment_components = []
            if len(control) != len(std_control) or len(std_control) != len(std_treatment):
                raise ValueError(f"{case}: standard-effect component counts do not match")
            for k, control_comp in enumerate(control):
                base = std_control[k]
                target = std_treatment[k]
                standard_treatment_components.append({
                    "component": control_comp["component"],
                    "weight": target["weight"],
                    "mu_R": control_comp["mu_R"] + target["mu_R"] - base["mu_R"],
                    "sigma_R": target["sigma_R"],
                    "mu_D": control_comp["mu_D"] + target["mu_D"] - base["mu_D"],
                    "sigma_D": target["sigma_D"],
                })
            for comp in control:
                all_parameters.append(output_row(
                    case, "H1-standard", "control", rho_label, rho, dict(comp),
                    "DERIVED_FROM_KM_FITTED_CONTROL", "KM-fitted H1-observed control",
                ))
            for comp in standard_treatment_components:
                all_parameters.append(output_row(
                    case, "H1-standard", "treatment", rho_label, rho, comp,
                    "DERIVED_LOCKED_STANDARD_EFFECT_FROM_KM_CONTROL", "E07 standard contrast + KM-fitted control",
                ))

    parameter_fields = list(all_parameters[0])
    with (OUT / "rho_sensitivity_km_target_parameters.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=parameter_fields)
        writer.writeheader()
        writer.writerows(all_parameters)
    metric_fields = list(metric_rows[0])
    with (OUT / "rho_sensitivity_km_target_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=metric_fields)
        writer.writeheader()
        writer.writerows(metric_rows)

    baseline_metrics = {
        (row["case"], row["arm"]): row
        for row in metric_rows
        if row["rho_label"] == "baseline"
    }
    acceptance_rows = []
    for row in metric_rows:
        baseline = baseline_metrics[(row["case"], row["arm"])]
        # KM targets contain digitization steps and are not expected to be
        # exactly generated by the restricted DGM family.  Use a conservative
        # absolute guard plus a relative-degradation guard versus the fitted
        # baseline rho.  The baseline itself is assessed by the same absolute
        # guard; low/high rho may not materially worsen that fit.
        primary_pass = (
            row["primary_weighted_rmse"] <= 0.03
            and row["primary_max_abs"] <= 0.08
            and row["primary_weighted_rmse"] <= baseline["primary_weighted_rmse"] + 0.01
            and row["primary_max_abs"] <= baseline["primary_max_abs"] + 0.02
        )
        os_pass = (
            row["death_check_weighted_rmse"] <= 0.03
            and row["death_check_max_abs"] <= 0.08
            and row["death_check_weighted_rmse"] <= baseline["death_check_weighted_rmse"] + 0.01
            and row["death_check_max_abs"] <= baseline["death_check_max_abs"] + 0.02
        )
        optimizer_pass = bool(row["optimizer_success"]) and bool(row["death_check_optimizer_success"])
        acceptance_rows.append({
            "case": row["case"], "scenario": row["scenario"], "arm": row["arm"],
            "rho_label": row["rho_label"], "rho": row["rho"],
            "primary_endpoint": row["primary_endpoint"],
            "primary_weighted_rmse": row["primary_weighted_rmse"],
            "primary_max_abs": row["primary_max_abs"],
            "baseline_primary_weighted_rmse": baseline["primary_weighted_rmse"],
            "baseline_primary_max_abs": baseline["primary_max_abs"],
            "os_weighted_rmse": row["os_weighted_rmse"],
            "os_max_abs": row["os_max_abs"],
            "death_check_weighted_rmse": row["death_check_weighted_rmse"],
            "death_check_max_abs": row["death_check_max_abs"],
            "baseline_death_check_weighted_rmse": baseline["death_check_weighted_rmse"],
            "baseline_death_check_max_abs": baseline["death_check_max_abs"],
            "pass_primary": primary_pass,
            "pass_os": os_pass,
            "pass_optimizer": optimizer_pass,
            "eligible_for_sensitivity_validation": primary_pass and os_pass and optimizer_pass,
            "reason": "PASS" if primary_pass and os_pass and optimizer_pass else "review fit or target mismatch",
        })
    with (OUT / "rho_sensitivity_km_target_acceptance.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(acceptance_rows[0]))
        writer.writeheader()
        writer.writerows(acceptance_rows)

    summary = {
        "status": "KM_TARGET_REFIT_COMPLETED_NO_MCMC",
        "target_manifest": str(TARGET_MANIFEST),
        "fit_definition": {
            "direct_fit": "H1-observed PFS/recurrence-CIF target plus OS check",
            "fixed": ["death marginal from E09 curve-fit parameters", "mixture weights", "model family", "trial process", "decision parameters"],
            "refit": ["mu_log_R", "sigma_log_R"],
            "derived": "H0 and H1-standard use KM-fitted control plus locked E07 standard treatment contrast",
        },
        "acceptance_rules": {
            "absolute_primary_weighted_rmse_max": 0.03,
            "absolute_primary_max_abs_max": 0.08,
            "absolute_os_weighted_rmse_max": 0.03,
            "absolute_os_max_abs_max": 0.08,
            "rho_vs_baseline_weighted_rmse_increase_max": 0.01,
            "rho_vs_baseline_max_abs_increase_max": 0.02,
        },
        "total_direct_groups": len(acceptance_rows),
        "eligible_direct_groups": sum(bool(r["eligible_for_sensitivity_validation"]) for r in acceptance_rows),
        "files": [
            "rho_sensitivity_km_target_parameters.csv",
            "rho_sensitivity_km_target_metrics.csv",
            "rho_sensitivity_km_target_acceptance.csv",
        ],
        "warning": "This output is deterministic DGM refitting only. It is not an MCMC independent-validation result and does not recalibrate decision thresholds.",
    }
    (OUT / "rho_sensitivity_km_target_fit_report.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"direct_groups": len(acceptance_rows), "eligible": summary["eligible_direct_groups"]}, ensure_ascii=False))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--case", nargs="+", choices=["survival_colon", "OAK", "IMpower131", "KEYNOTE-598"])
    parser.add_argument("--check-inputs", action="store_true", help="Check source inputs without optimization")
    args = parser.parse_args()
    TARGET_DIR = args.target_dir.resolve()
    TARGET_MANIFEST = TARGET_DIR / "target_manifest.json"
    SELECTED_CASES = set(args.case) if args.case else None
    if args.check_inputs:
        lock = json.loads(LOCK.read_text(encoding="utf-8"))
        targets = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))["target_rows"]
        rows = load_csv(DGM_INPUT)
        for spec in targets:
            path = Path(spec["target_file"])
            if not path.is_absolute(): path = TARGET_DIR / path
            if not path.is_file(): raise FileNotFoundError(path)
            curve = load_csv(path)
            if not curve: raise ValueError(f"Empty target: {path}")
        print(json.dumps({"status":"INPUTS_OK_NO_OPTIMIZATION", "targets":len(targets), "dgm_rows":len(rows), "cases":list(lock["rho_grid"])}))
    else:
        if args.output_dir is None: parser.error("--output-dir is required for fitting")
        OUT = require_external_output(args.output_dir, _SUPPORT_CODE_ROOT)
        if OUT == REPO_ROOT or REPO_ROOT in OUT.parents:
            parser.error("Choose an output directory outside the source package")
        OUT.mkdir(parents=True, exist_ok=False)
        main()
