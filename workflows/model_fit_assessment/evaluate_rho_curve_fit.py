#!/usr/bin/env python3
"""Evaluate the frozen primary and dependence-sensitivity DGMs.

All model parameters and observed inputs are read from the source package.
This entry evaluates curves; it does not fit parameters or simulate trials.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from functools import lru_cache
import math
from pathlib import Path
import sys
_SUPPORT_CODE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SUPPORT_CODE_ROOT / 'src/common'))
from output_paths import require_external_output

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr, ndtri


CASES = ("OAK", "IMpower131", "KEYNOTE-598", "survival_colon")
ARMS = ("control", "treatment")
MODELS = ("main_locked", "low", "high")
SCENARIO = "H1-observed"
DAYS_PER_MONTH = 365.25 / 12.0
RISK_WEIGHT_FLOOR = 0.025


@dataclass(frozen=True)
class Component:
    component: int
    weight: float
    mu_log_R: float
    sigma_log_R: float
    mu_log_D: float
    sigma_log_D: float
    rho: float
    cov_log_R_D: float


@dataclass
class EvaluationCurve:
    time: np.ndarray
    target: np.ndarray
    risk: np.ndarray
    support_end: float


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def require_file(root: Path, relative: str) -> Path:
    path = root / relative
    if not path.is_file():
        raise ValueError(f"Required package input is missing: {relative}")
    return path


def parameters(root: Path, selected_cases: list[str]) -> dict:
    formal = read_csv(require_file(root, "config/dgm_parameters.csv"))
    sensitivity = read_csv(require_file(
        root, "config/rho_parameters.csv"
    ))
    bundles = {}
    expected_components = {"OAK": 2, "KEYNOTE-598": 2,
                           "IMpower131": 1, "survival_colon": 1}
    for case in selected_cases:
        for model in MODELS:
            source = formal if model == "main_locked" else sensitivity
            for arm in ARMS:
                rows = [row for row in source
                        if row["case"] == case and row["scenario"] == SCENARIO
                        and row["arm"] == arm
                        and (model == "main_locked" or row["rho_label"] == model)]
                rows.sort(key=lambda row: int(row["component"]))
                if len(rows) != expected_components[case]:
                    raise ValueError(f"Unexpected component count: {case}/{model}/{arm}")
                components = []
                for row in rows:
                    component = Component(
                        component=int(row["component"]),
                        weight=float(row["weight"]),
                        mu_log_R=float(row["mu_log_R"]),
                        sigma_log_R=float(row["sigma_log_R"]),
                        mu_log_D=float(row["mu_log_D"]),
                        sigma_log_D=float(row["sigma_log_D"]),
                        rho=float(row["rho_log_R_D"] if model == "main_locked"
                                  else row["rho"]),
                        cov_log_R_D=float(row["cov_log_R_D"]),
                    )
                    values = list(asdict(component).values())
                    if not all(math.isfinite(value) for value in values):
                        raise ValueError(f"Non-finite parameter: {case}/{model}/{arm}")
                    if (component.weight <= 0 or component.sigma_log_R <= 0
                            or component.sigma_log_D <= 0 or abs(component.rho) >= 1):
                        raise ValueError(f"Invalid parameter: {case}/{model}/{arm}")
                    covariance = (component.rho * component.sigma_log_R
                                  * component.sigma_log_D)
                    if abs(covariance - component.cov_log_R_D) > 1e-8:
                        raise ValueError(f"Inconsistent covariance: {case}/{model}/{arm}")
                    components.append(component)
                if len({p.component for p in components}) != len(components):
                    raise ValueError(f"Duplicate components: {case}/{model}/{arm}")
                if abs(sum(p.weight for p in components) - 1.0) > 1e-8:
                    raise ValueError(f"Weights do not sum to one: {case}/{model}/{arm}")
                bundles[(case, model, arm)] = components
        for arm in ARMS:
            low = bundles[(case, "low", arm)]
            high = bundles[(case, "high", arm)]
            for first, second in zip(low, high):
                for name in ("weight", "mu_log_D", "sigma_log_D"):
                    if abs(getattr(first, name) - getattr(second, name)) > 1e-12:
                        raise ValueError(f"Death margin or weight changes across rho: {case}/{arm}")
    return bundles


def digitized_curve(path: Path) -> np.ndarray:
    values = np.loadtxt(path, delimiter=",", ndmin=2)
    if values.shape[1] != 2 or not np.all(np.isfinite(values)):
        raise ValueError(f"Invalid digitized curve: {path.name}")
    return values[np.argsort(values[:, 0], kind="stable")]


def lung_curves(root: Path, case: str) -> dict:
    prefix = "OAK_A_ITT" if case == "OAK" else (
        "KEYNOTE598" if case == "KEYNOTE-598" else case
    )
    input_dir = Path("inputs/clinical_data") / case
    curves = {}
    for endpoint in ("PFS", "OS"):
        for arm in ARMS:
            stem = f"{prefix}_{endpoint.lower()}_{arm}"
            coordinates = digitized_curve(require_file(root, str(input_dir / f"{stem}.csv")))
            suffix = "_pseudo_ipd.csv" if case == "KEYNOTE-598" else "_pseudo_ipd_local_20260731.csv"
            # The reconstructed OAK IPD uses the short OAK prefix.
            ipd_stem = stem.replace("OAK_A_ITT", "OAK")
            ipd = read_csv(require_file(root, str(Path("inputs/reconstructed_data/marginals") / case / f"{ipd_stem}{suffix}")))
            followup = np.sort(np.array([float(row["time_months"]) for row in ipd]))
            if len(followup) < 5 or not np.all(np.isfinite(followup)):
                raise ValueError(f"Invalid marginal IPD: {case}/{endpoint}/{arm}")
            if case == "OAK":
                end = 18.0 if endpoint == "PFS" else 27.0
            elif case == "IMpower131":
                end = 28.0 if endpoint == "PFS" else 36.0
            else:
                end = min(float(coordinates[-1, 0]), float(followup[-5]))
            grid = np.linspace(0.03, end - 0.0001, 1000)
            target = np.interp(grid, coordinates[:, 0], coordinates[:, 1])
            risk = len(followup) - np.searchsorted(followup, grid, side="left")
            if case == "KEYNOTE-598" and np.min(risk) < 5:
                raise ValueError(f"Unsupported KEYNOTE-598 evaluation point: {endpoint}/{arm}")
            curves[(arm, endpoint)] = EvaluationCurve(grid, target, risk, end)
    return curves


def observed_step_curve(times: np.ndarray, states: np.ndarray,
                        grid: np.ndarray, endpoint: str) -> tuple[np.ndarray, np.ndarray]:
    event_times = np.unique(times)
    values = []
    event_free = 1.0
    value = 0.0 if endpoint == "recurrence_CIF" else 1.0
    for time in event_times:
        risk = np.count_nonzero(times >= time)
        if endpoint == "recurrence_CIF":
            recurrence = np.count_nonzero((times == time) & (states == 1))
            death = np.count_nonzero((times == time) & (states == 2))
            value += event_free * recurrence / risk
            event_free *= 1.0 - (recurrence + death) / risk
        else:
            death = np.count_nonzero((times == time) & (states == 1))
            value *= 1.0 - death / risk
        values.append(value)
    index = np.searchsorted(event_times, grid, side="right") - 1
    initial = 0.0 if endpoint == "recurrence_CIF" else 1.0
    target = np.where(index >= 0, np.array(values)[np.maximum(index, 0)], initial)
    risk = len(times) - np.searchsorted(np.sort(times), grid, side="left")
    return target, risk


def colon_curves(root: Path) -> dict:
    path = require_file(root, "inputs/clinical_data/survival_colon.csv")
    patients = {}
    for row in read_csv(path):
        if row["rx"] not in ("Obs", "Lev+5FU"):
            continue
        endpoint = int(row["etype"])
        patient = patients.setdefault(int(row["id"]), {})
        if endpoint in patient:
            raise ValueError("Duplicate patient endpoint in survival_colon input")
        patient[endpoint] = row
    if len(patients) != 619:
        raise ValueError("Expected 619 independent survival_colon patients")
    records = {arm: [] for arm in ARMS}
    for patient in patients.values():
        if set(patient) != {1, 2} or patient[1]["rx"] != patient[2]["rx"]:
            raise ValueError("Incomplete or inconsistent survival_colon patient endpoints")
        recurrence, death = patient[1], patient[2]
        arm = "control" if recurrence["rx"] == "Obs" else "treatment"
        records[arm].append((float(recurrence["time"]) / DAYS_PER_MONTH,
                             float(death["time"]) / DAYS_PER_MONTH,
                             int(recurrence["status"]), int(death["status"])))
    grid = np.arange(0.0, 84.0001, 0.5)
    curves = {}
    for arm in ARMS:
        if len(records[arm]) != (315 if arm == "control" else 304):
            raise ValueError(f"Unexpected survival_colon arm size: {arm}")
        recurrence, death, recurrence_status, death_status = np.array(records[arm]).T
        first = np.minimum(recurrence, death)
        # Same-day death takes precedence, following the original Figure S1 rule.
        state = np.where((recurrence_status == 1) & (recurrence < death), 1,
                         np.where((death_status == 1) & (death <= recurrence), 2, 0))
        for endpoint, time, status in (("recurrence_CIF", first, state),
                                       ("OS", death, death_status)):
            target, risk = observed_step_curve(time, status, grid, endpoint)
            curves[(arm, endpoint)] = EvaluationCurve(grid, target, risk, 84.0)
    return curves


@lru_cache(maxsize=None)
def quadrature(nodes: int) -> tuple[np.ndarray, np.ndarray]:
    return np.polynomial.legendre.leggauss(nodes)


def pfs_component(time: np.ndarray, parameter: Component, nodes: int) -> np.ndarray:
    result = np.ones_like(time, dtype=float)
    positive = time > 0
    log_time = np.log(time[positive])
    a = (log_time - parameter.mu_log_R) / parameter.sigma_log_R
    b = (log_time - parameter.mu_log_D) / parameter.sigma_log_D
    x, weight = quadrature(nodes)
    half_width = ndtr(-b) / 2.0
    uniform = 1.0 - half_width[:, None] + half_width[:, None] * x
    normal = ndtri(np.clip(uniform, 1e-14, 1.0 - 1e-14))
    result[positive] = half_width * (
        ndtr((parameter.rho * normal - a[:, None])
             / math.sqrt(1.0 - parameter.rho**2)) @ weight
    )
    return np.clip(result, 0.0, 1.0)


def recurrence_component(time: np.ndarray, parameter: Component) -> np.ndarray:
    """Pr(T_R <= t, T_R < T_D), with death as a competing event."""
    result = np.zeros_like(time, dtype=float)
    denominator = parameter.sigma_log_D * math.sqrt(1.0 - parameter.rho**2)
    for index, value in enumerate(time):
        if value <= 0:
            continue
        bound = (math.log(value) - parameter.mu_log_R) / parameter.sigma_log_R
        integral, error = quad(
            lambda z: math.exp(-0.5 * z*z) / math.sqrt(2.0 * math.pi)
            * ndtr((parameter.mu_log_D - parameter.mu_log_R
                    + (parameter.rho * parameter.sigma_log_D
                       - parameter.sigma_log_R) * z) / denominator),
            -np.inf, bound, epsabs=2e-12, epsrel=2e-11, limit=200,
        )
        if error > 1e-8:
            raise ValueError("Recurrence CIF integration did not reach the requested precision")
        result[index] = integral
    return np.clip(result, 0.0, 1.0)


def theoretical_curve(time: np.ndarray, components: list[Component],
                      endpoint: str, nodes: int) -> np.ndarray:
    result = np.zeros_like(time, dtype=float)
    for parameter in components:
        if endpoint == "PFS":
            value = pfs_component(time, parameter, nodes)
        elif endpoint == "recurrence_CIF":
            value = recurrence_component(time, parameter)
        else:
            value = np.ones_like(time, dtype=float)
            positive = time > 0
            value[positive] = ndtr((parameter.mu_log_D - np.log(time[positive]))
                                   / parameter.sigma_log_D)
        result += parameter.weight * value
    return result


def error_metrics(curve: EvaluationCurve, fitted: np.ndarray) -> dict:
    delta = fitted - curve.target
    weights = np.maximum(curve.risk / max(curve.risk[0], 1), RISK_WEIGHT_FLOOR)
    maximum = int(np.argmax(np.abs(delta)))
    return dict(mae=float(np.mean(np.abs(delta))),
                rmse=float(np.sqrt(np.mean(delta**2))),
                risk_weighted_rmse=float(np.sqrt(np.average(delta**2, weights=weights))),
                max_abs_error=float(np.max(np.abs(delta))),
                time_at_max_error_months=float(curve.time[maximum]))


def reference_comparisons(root: Path, curves: dict, fits: dict) -> list[dict]:
    """Compare primary curves with supplied fixed curve references."""
    rows = []
    sources = [
        (root / "inputs/clinical_data/evaluation_curve_sources.csv", None),
    ]
    for path, fixed_case in sources:
        if not path.is_file():
            continue
        reference = read_csv(path)
        for key, curve in curves.items():
            case, arm, endpoint = key
            if fixed_case is not None and case != fixed_case:
                continue
            source = [r for r in reference
                      if (fixed_case is not None or r["case"] == case)
                      and r["arm"] == arm and r["endpoint"] == endpoint
                      and ("layer" not in r or r["layer"] == "final_dgm_fit")]
            if not source:
                continue
            source.sort(key=lambda r: float(r["time_months"]))
            time = np.array([float(r["time_months"]) for r in source])
            target = np.array([float(r["target"]) for r in source])
            risk = np.array([float(r["risk"]) for r in source])
            if (len(time) != len(curve.time) or not np.allclose(time, curve.time, rtol=0, atol=1e-12)
                    or not np.allclose(target, curve.target, rtol=0, atol=2e-14)
                    or not np.array_equal(risk, curve.risk)):
                raise ValueError(f"Fixed curve reference differs from observed inputs: {key}")
            original = np.array([float(r["fitted"]) for r in source])
            primary = fits[(case, "main_locked", arm, endpoint)]
            old_metrics = error_metrics(curve, original)
            new_metrics = error_metrics(curve, primary)
            row = dict(case=case, arm=arm, endpoint=endpoint,
                       reference_type="mean_of_1000_generated_curves" if fixed_case else "locked_figure_curve",
                       max_abs_primary_vs_reference=float(np.max(np.abs(primary - original))))
            for name in ("mae", "rmse", "risk_weighted_rmse", "max_abs_error"):
                row[f"primary_{name}"] = new_metrics[name]
                row[f"reference_{name}"] = old_metrics[name]
            rows.append(row)
    return rows


def run(package_root: Path, output_dir: Path, selected_cases: list[str],
        nodes: int, write_curves: bool) -> None:
    bundles = parameters(package_root, selected_cases)
    metrics, windows, used, points = [], [], [], []
    observed, fitted = {}, {}
    for case in selected_cases:
        curves = colon_curves(package_root) if case == "survival_colon" else lung_curves(package_root, case)
        for (arm, endpoint), curve in curves.items():
            observed[(case, arm, endpoint)] = curve
            windows.append(dict(case=case, arm=arm, endpoint=endpoint,
                                start_months=float(curve.time[0]),
                                support_end_months=curve.support_end,
                                last_grid_months=float(curve.time[-1]),
                                n_grid=len(curve.time), minimum_risk=int(np.min(curve.risk)),
                                risk_at_start=int(curve.risk[0])))
            for model in MODELS:
                components = bundles[(case, model, arm)]
                value = theoretical_curve(curve.time, components, endpoint, nodes)
                fitted[(case, model, arm, endpoint)] = value
                metrics.append(dict(case=case, model=model, scenario=SCENARIO,
                                    effective_rho=components[0].rho, arm=arm, endpoint=endpoint,
                                    start_months=float(curve.time[0]),
                                    last_months=float(curve.time[-1]), n_grid=len(curve.time),
                                    **error_metrics(curve, value)))
                if write_curves:
                    for time, target, fit, risk in zip(curve.time, curve.target, value, curve.risk):
                        points.append(dict(case=case, model=model, scenario=SCENARIO,
                                           effective_rho=components[0].rho, arm=arm, endpoint=endpoint,
                                           time_months=float(time), target=float(target),
                                           fitted=float(fit), risk=int(risk)))
        for model in MODELS:
            for arm in ARMS:
                used.extend(dict(case=case, model=model, scenario=SCENARIO, arm=arm,
                                 **asdict(component)) for component in bundles[(case, model, arm)])
    references = reference_comparisons(package_root, observed, fitted)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "curve_fit_metrics.csv", metrics)
    write_csv(output_dir / "curve_fit_windows.csv", windows)
    write_csv(output_dir / "used_parameters.csv", used)
    write_csv(output_dir / "primary_reference_comparison.csv", references)
    if write_curves:
        write_csv(output_dir / "evaluation_curves.csv", points)
    print(f"Evaluated {len(metrics)} curves. Results: {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path,
                        default=next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir()),
                        help="Source-package root; defaults to this script's package.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Directory outside the source package for generated results.")
    parser.add_argument("--case", choices=("all",) + CASES, default="all")
    parser.add_argument("--quadrature-nodes", type=int, default=640,
                        help="Gauss-Legendre nodes for lung-cancer PFS curves (default: 640).")
    parser.add_argument("--write-curves", action="store_true",
                        help="Also save the evaluation-point targets, fitted values and risk counts.")
    args = parser.parse_args()
    root = args.package_root.expanduser().resolve()
    output = require_external_output(args.output_dir, _SUPPORT_CODE_ROOT)
    if output == root or root in output.parents:
        parser.error("--output-dir must be outside the source package")
    if not root.is_dir():
        parser.error("--package-root must be an existing source-package directory")
    if args.quadrature_nodes < 32:
        parser.error("--quadrature-nodes must be at least 32")
    selected = list(CASES) if args.case == "all" else [args.case]
    try:
        run(root, output, selected, args.quadrature_nodes, args.write_curves)
    except (ValueError, OSError, KeyError) as error:
        parser.exit(1, f"Input or calculation error: {error}\n")


if __name__ == "__main__":
    main()
