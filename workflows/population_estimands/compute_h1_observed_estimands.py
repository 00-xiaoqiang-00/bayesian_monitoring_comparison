"""Evaluate the locked H1-observed population estimands.

The numerical functions are unchanged from the source identified in the saved
computation configuration. Default execution verifies the existing archive;
only explicit --recompute performs numerical integration / Monte Carlo work.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
import sys
_SOURCE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SOURCE_ROOT / 'src/dgm'))
from time import perf_counter

import statistics
import sys

import compute_dgm_estimands as dgm

PACKAGE = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
DEFAULT_ARCHIVE = PACKAGE.parent / "04_Data_and_Results/Population_Estimands/Final/H1_Observed"


def load_numerical_runtime():
    global np, scipy, quad, brentq, multivariate_normal, norm, student_t
    import numpy as np
    import scipy
    from scipy.integrate import quad
    from scipy.optimize import brentq
    from scipy.stats import multivariate_normal, norm, t as student_t
    dgm.np = np
    dgm.brentq = brentq
    dgm.multivariate_normal = multivariate_normal
    dgm.norm = norm

WR_BATCHES = 100
WR_PAIRS = 200_000
COX_BATCHES = 24
COX_N_ARM = 250_000
Z_BOUND = 12.0
QUAD_ABS = 1e-9
QUAD_REL = 1e-10


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def verify_inputs(input_dir: Path, specs: list) -> dict:
    parameters_path = input_dir / "dgm_parameters.csv"
    trials_path = input_dir / "trial_process.csv"
    actual_parameters = read_csv(parameters_path)
    expected_parameters = []
    for scenario in ("standard", "observed"):
        expected_parameters.extend(dgm.parameter_rows(specs, scenario))
    for row in dgm.parameter_rows(specs, "standard"):
        if row["arm"] == "control":
            expected_parameters.append(dict(row, scenario="H0", arm="control"))
            expected_parameters.append(dict(row, scenario="H0", arm="treatment"))

    def parameter_key(row):
        return row["case"], row["scenario"], row["arm"], int(row["component"])

    actual_map = {parameter_key(row): row for row in actual_parameters}
    expected_map = {parameter_key(row): row for row in expected_parameters}
    if len(actual_map) != len(actual_parameters) or set(actual_map) != set(expected_map):
        raise ValueError("DGM input keys differ from the locked source, or are duplicated.")
    parameter_fields = 0
    for key, expected in expected_map.items():
        actual = actual_map[key]
        if set(actual) != set(expected):
            raise ValueError(f"DGM columns differ for {key}.")
        for field, value in expected.items():
            if isinstance(value, (int, float)):
                matched = math.isclose(float(actual[field]), value, rel_tol=1e-9, abs_tol=1e-9)
            else:
                matched = actual[field] == str(value)
            if not matched:
                raise ValueError(f"DGM mismatch: {key}, {field}.")
            parameter_fields += 1
    actual_trials = read_csv(trials_path)
    trials_map = {row["case"]: row for row in actual_trials}
    expected_trials = dgm.trial_rows(specs)
    if len(trials_map) != len(actual_trials) or set(trials_map) != {row["case"] for row in expected_trials}:
        raise ValueError("Trial-process cases differ from the locked source, or are duplicated.")
    trial_fields = 0
    for expected in expected_trials:
        actual = trials_map[expected["case"]]
        if set(actual) != set(expected):
            raise ValueError(f"Trial-process columns differ for {expected['case']}.")
        for field, value in expected.items():
            if isinstance(value, (int, float)):
                matched = math.isclose(float(actual[field]), value, rel_tol=1e-12, abs_tol=1e-12)
            else:
                matched = actual[field] == str(value)
            if not matched:
                raise ValueError(f"Trial-process mismatch: {expected['case']}, {field}.")
            trial_fields += 1
    return {
        "status": "passed",
        "DGM_rows": len(actual_parameters),
        "DGM_fields_matched": parameter_fields,
        "trial_rows": len(actual_trials),
        "trial_fields_matched": trial_fields,
        "input_sha256": {path.name: sha256(path) for path in (parameters_path, trials_path)},
        "precision_note": "Numerical values use the locked Python specification after checking the exported CSVs. The trial CSV's censoring description is formatted; KEYNOTE-598 uses the full stored bound 141.8559446075.",
    }


def lognormal_restricted_mean(mu: float, sigma: float, horizon: float) -> float:
    if horizon <= 0:
        return 0.0
    log_horizon = math.log(horizon)
    return float(
        math.exp(mu + 0.5 * sigma * sigma)
        * norm.cdf((log_horizon - mu - sigma * sigma) / sigma)
        + horizon * norm.cdf((mu - log_horizon) / sigma)
    )


def checked_quad(function, lower: float, upper: float, *, epsabs=QUAD_ABS, epsrel=QUAD_REL):
    if lower >= upper:
        return 0.0, 0.0
    answer = quad(function, lower, upper, epsabs=epsabs, epsrel=epsrel, limit=250, full_output=1)
    if len(answer) != 3:
        raise ArithmeticError("Numerical quadrature did not converge: " + str(answer[3]))
    value, error, _ = answer
    if not (math.isfinite(value) and math.isfinite(error)):
        raise ArithmeticError("Non-finite quadrature result.")
    return float(value), float(error)


def component_composite_mean(component, horizon: float) -> tuple[float, dict]:
    conditional_sigma = component.sigma_d * math.sqrt(1.0 - component.rho ** 2)

    def integrand(z):
        recurrence = math.exp(component.mu_r + component.sigma_r * z)
        conditional_mu = component.mu_d + component.rho * component.sigma_d * z
        return float(norm.pdf(z)) * lognormal_restricted_mean(
            conditional_mu, conditional_sigma, min(recurrence, horizon)
        )

    split = (math.log(horizon) - component.mu_r) / component.sigma_r
    boundaries = [-Z_BOUND]
    if -Z_BOUND < split < Z_BOUND:
        boundaries.append(split)
    boundaries.append(Z_BOUND)
    total = 0.0
    error_total = 0.0
    for lower, upper in zip(boundaries[:-1], boundaries[1:]):
        value, error = checked_quad(integrand, lower, upper)
        total += value
        error_total += error
    return total, {
        "estimated_quadrature_absolute_error_months": error_total,
        "omitted_tail_upper_bound_months": float(2.0 * horizon * norm.sf(Z_BOUND)),
        "z_split_at_R_equals_tau": split,
        "integration_segments": boundaries,
    }


def arm_composite_mean(arm, horizon: float) -> tuple[float, dict]:
    values = []
    value_total = 0.0
    error_total = 0.0
    tail_total = 0.0
    for index, component in enumerate(arm.components, 1):
        value, detail = component_composite_mean(component, horizon)
        value_total += component.weight * value
        error_total += component.weight * detail["estimated_quadrature_absolute_error_months"]
        tail_total += component.weight * detail["omitted_tail_upper_bound_months"]
        values.append(dict(component=index, weight=component.weight, restricted_mean_months=value, **detail))
    return value_total, {
        "estimated_quadrature_absolute_error_months": error_total,
        "omitted_tail_upper_bound_months": tail_total,
        "components": values,
    }


def arm_death_mean(arm, horizon: float) -> float:
    return sum(c.weight * lognormal_restricted_mean(c.mu_d, c.sigma_d, horizon) for c in arm.components)


def crosscheck_component_survival_integral(component, horizon: float, reference: float) -> dict:
    # Independent order of integration: integrate P(R>t,D>t) over t.
    conditional_sigma = component.sigma_d * math.sqrt(1.0 - component.rho ** 2)
    maximum_inner_error = 0.0

    def survival(time):
        nonlocal maximum_inner_error
        if time <= 0:
            return 1.0
        log_time = math.log(time)
        lower = max(-Z_BOUND, (log_time - component.mu_r) / component.sigma_r)
        if lower >= Z_BOUND:
            return 0.0

        def integrand(z):
            conditional_mu = component.mu_d + component.rho * component.sigma_d * z
            return float(norm.pdf(z) * norm.sf((log_time - conditional_mu) / conditional_sigma))

        value, error = checked_quad(integrand, lower, Z_BOUND, epsabs=1e-10, epsrel=1e-10)
        maximum_inner_error = max(maximum_inner_error, error)
        return value

    value, outer_error = checked_quad(survival, 0.0, horizon, epsabs=1e-7, epsrel=1e-8)
    difference = abs(value - reference)
    tolerance = max(2e-5, 20.0 * (outer_error + horizon * maximum_inner_error))
    if difference > tolerance:
        raise ArithmeticError(f"Independent survival-integral check failed: difference={difference}, tolerance={tolerance}.")
    return {
        "method": "Integral from 0 to tau of conditional-normal composite survival; control component 1",
        "restricted_mean_months": value,
        "conditional_expectation_reference_months": reference,
        "absolute_difference_months": difference,
        "outer_quadrature_error_estimate_months": outer_error,
        "maximum_inner_probability_error_estimate": maximum_inner_error,
        "acceptance_tolerance_months": tolerance,
        "status": "passed",
    }


def deterministic_estimands(spec) -> tuple[dict, dict]:
    horizon_rmst = spec.trial.tau_rmst
    horizon_rmut = spec.trial.tau_rmut
    control_rmst, control_rmst_detail = arm_composite_mean(spec.control, horizon_rmst)
    observed_rmst, observed_rmst_detail = arm_composite_mean(spec.observed, horizon_rmst)
    if horizon_rmut == horizon_rmst:
        control_composite, control_composite_detail = control_rmst, control_rmst_detail
        observed_composite, observed_composite_detail = observed_rmst, observed_rmst_detail
    else:
        control_composite, control_composite_detail = arm_composite_mean(spec.control, horizon_rmut)
        observed_composite, observed_composite_detail = arm_composite_mean(spec.observed, horizon_rmut)
    control_death = arm_death_mean(spec.control, horizon_rmut)
    observed_death = arm_death_mean(spec.observed, horizon_rmut)
    delta_rmst = observed_rmst - control_rmst
    delta_sc = observed_composite - control_composite
    delta_sd = observed_death - control_death
    weight = spec.trial.w_r
    delta_rmut = (1.0 - weight) * delta_sd + weight * delta_sc
    rmst_error = sum(x["estimated_quadrature_absolute_error_months"] for x in (control_rmst_detail, observed_rmst_detail))
    rmst_tail = sum(x["omitted_tail_upper_bound_months"] for x in (control_rmst_detail, observed_rmst_detail))
    rmut_error = weight * sum(x["estimated_quadrature_absolute_error_months"] for x in (control_composite_detail, observed_composite_detail))
    rmut_tail = weight * sum(x["omitted_tail_upper_bound_months"] for x in (control_composite_detail, observed_composite_detail))
    first_reference = control_rmst_detail["components"][0]["restricted_mean_months"]
    crosscheck = crosscheck_component_survival_integral(spec.control.components[0], horizon_rmst, first_reference)
    return {
        "delta_rmst_months": delta_rmst,
        "tau_rmst_months": horizon_rmst,
        "delta_rmst_quadrature_error_estimate_months": rmst_error,
        "delta_rmst_omitted_tail_bound_months": rmst_tail,
        "delta_rmut_months": delta_rmut,
        "tau_rmut_months": horizon_rmut,
        "delta_rmut_quadrature_error_estimate_months": rmut_error,
        "delta_rmut_omitted_tail_bound_months": rmut_tail,
        "rmut_weight_recurrence_state": weight,
        "delta_integral_sd_months": delta_sd,
        "delta_integral_sc_months": delta_sc,
        "control_composite_rmst_months": control_rmst,
        "observed_composite_rmst_months": observed_rmst,
        "control_composite_integral_at_rmut_tau_months": control_composite,
        "observed_composite_integral_at_rmut_tau_months": observed_composite,
        "control_death_integral_at_rmut_tau_months": control_death,
        "observed_death_integral_at_rmut_tau_months": observed_death,
    }, {
        "control_RMST": control_rmst_detail,
        "observed_RMST": observed_rmst_detail,
        "control_composite_at_RMUT_tau": control_composite_detail,
        "observed_composite_at_RMUT_tau": observed_composite_detail,
        "independent_survival_integral_crosscheck": crosscheck,
        "error_note": "Quadrature errors are numerical estimates. Omitted-normal-tail bounds are analytic. Floating-point rounding of closed-form lognormal restricted means is not included in these estimates.",
    }


def wr_conditional_batch(spec, rng) -> dict:
    rt, dt = dgm.draw_arm(spec.observed, WR_PAIRS, rng)
    rc, dc = dgm.draw_arm(spec.control, WR_PAIRS, rng)
    death_min = np.minimum(dt, dc)
    death_win = (dc < dt) & (dc < spec.trial.final_analysis_time)
    death_loss = (dt < dc) & (dt < spec.trial.final_analysis_time)
    recurrence_win = (rc < np.minimum(rt, death_min)) & (rc < spec.trial.final_analysis_time)
    recurrence_loss = (rt < np.minimum(rc, death_min)) & (rt < spec.trial.final_analysis_time)
    rec_win_contribution = np.zeros(WR_PAIRS, dtype=float)
    rec_loss_contribution = np.zeros(WR_PAIRS, dtype=float)
    for mask, recurrence, contribution in (
        (recurrence_win, rc, rec_win_contribution),
        (recurrence_loss, rt, rec_loss_contribution),
    ):
        at_recurrence = dgm.censor_survival(recurrence[mask], spec.trial)
        at_death = dgm.censor_survival(death_min[mask], spec.trial)
        if np.any(at_recurrence <= 0) or np.any(at_death > at_recurrence + 1e-12):
            raise ArithmeticError("Invalid censoring-survival ratio in conditional WR.")
        values = 1.0 - np.square(at_death / at_recurrence)
        if np.any(values < -1e-12) or np.any(values > 1 + 1e-12):
            raise ArithmeticError("Conditional WR recurrence contribution is outside [0,1].")
        contribution[mask] = np.clip(values, 0.0, 1.0)
    wins = death_win.astype(float) + rec_win_contribution
    losses = death_loss.astype(float) + rec_loss_contribution
    return {
        "pairs": WR_PAIRS,
        "conditional_weighted_wins_mean": float(np.mean(wins)),
        "conditional_weighted_losses_mean": float(np.mean(losses)),
        "expected_ipcw_death_win_per_pair": float(np.mean(death_win)),
        "expected_ipcw_death_loss_per_pair": float(np.mean(death_loss)),
        "expected_ipcw_recurrence_win_per_pair": float(np.mean(rec_win_contribution)),
        "expected_ipcw_recurrence_loss_per_pair": float(np.mean(rec_loss_contribution)),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_case(spec, index: int, output: Path) -> tuple[dict, dict, list, list]:
    started = perf_counter()
    seed = 2026102300 + index * 100
    print(f"{spec.case}: deterministic RMST/RMUT integration", flush=True)
    integral_result, integral_detail = deterministic_estimands(spec)
    observed_spec = replace(spec, standard=spec.observed)
    cox_rows = []
    wr_rows = []
    for batch in range(COX_BATCHES):
        rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(0, batch)))
        hr = dgm.cox_hr_one_batch(observed_spec, COX_N_ARM, rng)
        if not math.isfinite(hr) or hr <= 0:
            raise ArithmeticError(f"Invalid Cox HR: {spec.case}, batch {batch + 1}.")
        cox_rows.append({"case": spec.case, "scenario": "H1-observed", "batch": batch + 1,
                         "seed": seed, "spawn_key": f"0;{batch}", "n_per_arm": COX_N_ARM,
                         "cox_hr": hr, "cox_log_hr": math.log(hr)})
        if (batch + 1) % 6 == 0:
            print(f"{spec.case}: Cox {batch + 1}/{COX_BATCHES} batches", flush=True)
    for batch in range(WR_BATCHES):
        rng = np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(1, batch)))
        wr_rows.append({"case": spec.case, "scenario": "H1-observed", "batch": batch + 1,
                        "seed": seed, "spawn_key": f"1;{batch}", **wr_conditional_batch(spec, rng)})
        if (batch + 1) % 20 == 0:
            print(f"{spec.case}: WR {batch + 1}/{WR_BATCHES} batches", flush=True)

    log_hrs = np.array([row["cox_log_hr"] for row in cox_rows])
    log_hr = float(np.mean(log_hrs))
    log_hr_mcse = float(np.std(log_hrs, ddof=1) / math.sqrt(COX_BATCHES))
    hr = math.exp(log_hr)
    hr_critical = float(student_t.ppf(0.975, COX_BATCHES - 1))
    wins = np.array([row["conditional_weighted_wins_mean"] for row in wr_rows])
    losses = np.array([row["conditional_weighted_losses_mean"] for row in wr_rows])
    wr = float(np.mean(wins) / np.mean(losses))
    wr_mcse = float(np.std(wins - wr * losses, ddof=1) / math.sqrt(WR_BATCHES) / np.mean(losses))
    wr_critical = float(student_t.ppf(0.975, WR_BATCHES - 1))
    result = {
        "case": spec.case,
        "scenario": "H1-observed",
        "locked_dgm_version_identifier": spec.version,
        "generator": spec.generator,
        "wr": wr,
        "wr_mcse": wr_mcse,
        "wr_numerical_95_lower": wr - wr_critical * wr_mcse,
        "wr_numerical_95_upper": wr + wr_critical * wr_mcse,
        "wr_log": math.log(wr),
        "wr_expected_ipcw_wins_per_pair": float(np.mean(wins)),
        "wr_expected_ipcw_losses_per_pair": float(np.mean(losses)),
        "cox_hr": hr,
        "cox_hr_mcse": hr * log_hr_mcse,
        "cox_log_hr": log_hr,
        "cox_log_hr_mcse": log_hr_mcse,
        "cox_hr_numerical_95_lower": math.exp(log_hr - hr_critical * log_hr_mcse),
        "cox_hr_numerical_95_upper": math.exp(log_hr + hr_critical * log_hr_mcse),
        **integral_result,
        "cox_and_wr_final_trial_process": (
            f"Nmax={spec.trial.nmax}; recruitment={spec.trial.recruitment_duration}; "
            f"final={spec.trial.final_analysis_time}; random_censor=Uniform(0,{spec.trial.cens_upper:.13g})"
        ),
        "cox_batches": COX_BATCHES,
        "cox_n_per_arm_per_batch": COX_N_ARM,
        "wr_batches": WR_BATCHES,
        "wr_pairs_per_batch": WR_PAIRS,
        "seed": seed,
    }
    print(f"{spec.case}: completed in {perf_counter() - started:.1f} s; WR={wr:.8f} (MCSE {wr_mcse:.3g}), Cox HR={hr:.8f} (MCSE {hr * log_hr_mcse:.3g})", flush=True)
    return result, integral_detail, cox_rows, wr_rows


def historical_recompute(argv=None) -> None:
    parser = argparse.ArgumentParser(description="Compute locked H1-observed DGM estimands on the analysis server.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cases", nargs="+", choices=[s.case for s in dgm.build_specs()])
    args = parser.parse_args(argv)
    input_dir = args.input_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        parser.error("The output directory must not already exist.")
    if args.cases is not None and len(args.cases) != len(set(args.cases)):
        parser.error("Do not repeat case names.")
    specs = dgm.build_specs()
    verified = verify_inputs(input_dir, specs)
    output.mkdir(parents=True)
    configuration = {
        "scenario": "H1-observed",
        "started_utc": utc_now(),
        "status": "running",
        "input_verification": verified,
        "source_sha256": {Path(__file__).name: sha256(Path(__file__)), Path(dgm.__file__).name: sha256(Path(dgm.__file__))},
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform()},
        "runtime_note": "Runtime information records the actual server environment. It is provenance for this numerical computation, not a change to the manuscript's R 4.4.x environment statement.",
        "definitions": {
            "WR": "Death-prioritized IPCW population pairwise win ratio under the locked final trial process. The independent censoring variables are integrated out analytically conditional on latent event times; this preserves the original target.",
            "Cox_HR": "Large-sample approximation to the equal-allocation marginal Cox partial-likelihood projection for min(R,D), under the locked final trial process. This need not be a time-constant true hazard ratio.",
            "RMST": "Observed-treatment minus control integral of P(min(R,D)>t), from zero through the method-specific tau.",
            "RMUT": "Observed-treatment minus control integral of (1-w_R)P(D>t)+w_R P(min(R,D)>t), with w_R=0.5 and the locked RMUT tau.",
            "numerical_intervals": "Monte Carlo numerical uncertainty only; no original-trial sampling, DGM-fit, or structural-assumption uncertainty is included.",
        },
        "numerical_configuration": {
            "base_seed": "2026102300 + 100 * index in locked case order",
            "independent_rng_streams": "NumPy SeedSequence(base_seed, spawn_key=(method_index, zero_based_batch)); Cox method_index=0, WR=1",
            "cox_batches": COX_BATCHES, "cox_n_per_arm_per_batch": COX_N_ARM,
            "wr_batches": WR_BATCHES, "wr_pairs_per_batch": WR_PAIRS,
            "wr_mcse": "sd(W_batch_mean - WR*L_batch_mean) / sqrt(B) / mean(L_batch_mean); ratio delta method",
            "wr_interval": "WR +/- t_0.975,99 * MCSE",
            "cox_point_estimate": "exp(mean(log(HR_batch)))",
            "cox_interval": "exp(mean(log(HR_batch)) +/- t_0.975,23 * sd(log(HR_batch))/sqrt(24))",
            "cox_precision_limitation": "Batch Monte Carlo SE does not measure the finite-n approximation bias of Cox fits; comparison with larger per-batch sample sizes is a separate numerical check.",
            "normal_integration_z_bound": Z_BOUND, "quadrature_epsabs": QUAD_ABS, "quadrature_epsrel": QUAD_REL,
            "closed_form_death_restricted_mean": True,
            "wr_no_weight_truncation": True,
        },
        "implementation_note": "This computes H1-observed using a new numerical implementation and new independent random streams. It does not claim bitwise reproduction of the archived H1-standard calculation. The locked event distributions and estimand definitions are unchanged.",
        "locked_version_identifier_note": "Version/generator strings identify the pre-existing locked specification. In particular, OAK's equal_H1_variance label concerns H1-standard; H1-observed uses its own stored unequal dispersions.",
        "case_configurations": [],
    }
    config_path = output / "computation_config.json"

    def save_config():
        config_path.write_text(json.dumps(configuration, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    save_config()
    results = []
    integral_diagnostics = {}
    all_cox_batches = []
    all_wr_batches = []
    try:
        for index, spec in enumerate(specs):
            if args.cases is not None and spec.case not in args.cases:
                continue
            row, diagnostics, cox_rows, wr_rows = run_case(spec, index, output)
            results.append(row)
            integral_diagnostics[spec.case] = diagnostics
            all_cox_batches.extend(cox_rows)
            all_wr_batches.extend(wr_rows)
            configuration["case_configurations"].append({
                "case": spec.case, "seed": row["seed"],
                "trial": dgm.trial_rows([spec])[0],
                "H1_observed_parameters": dgm.parameter_rows([spec], "observed"),
            })
            write_csv(output / "H1_observed_true_estimands.csv", results)
            save_config()
        write_csv(output / "cox_batches.csv", all_cox_batches)
        write_csv(output / "wr_batches.csv", all_wr_batches)
        (output / "integration_diagnostics.json").write_text(json.dumps(integral_diagnostics, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        readme = """# H1-observed population estimands

These results evaluate the locked H1-observed experimental distribution against its locked control distribution. They describe the simulation model, rather than an unknown clinical population effect. No monitoring-rule calibration or operating-characteristic simulation was rerun.

`H1_observed_true_estimands.csv` contains the estimands and numerical uncertainty. Values have computational precision, not corresponding clinical accuracy. WR and Cox intervals describe Monte Carlo uncertainty only. RMST/RMUT quadrature errors are estimates; the omitted-normal-tail bounds are analytic. The original-data sampling and fitted-model uncertainty are not included.

WR preserves the original death-prioritized IPCW target and final follow-up process. Conditional expectation integrates out independent censoring analytically. For mD=min(Dt,Dc), a death win contributes I(Dc<Dt,Dc<T); a recurrence win contributes I(Rc<min(Rt,Dt,Dc),Rc<T) * [1-(G(mD)/G(Rc))^2]. Losses are symmetric. Conditional contributions may overlap for a latent pair because they average over different censoring realizations; they are not raw mutually exclusive outcome probabilities. No weights or follow-up horizons were truncated.

Cox uses the first event min(R,D), equal allocation, staggered recruitment, random independent censoring and the locked final calendar cutoff. Under nonproportional hazards this is a process-dependent projection. It is not the published OS HR. Cox uses the locked fitting function with H1-observed substituted for the treatment arm. The result is a large-sample projection approximation; its batch Monte Carlo SE does not include finite-n approximation bias.

RMST integrates composite-event survival. RMUT uses (1-w_R) times death survival plus w_R times composite-event survival, with w_R=0.5. Their units are months. The KEYNOTE-598 time windows differ: RMST 12 months, RMUT 18 months; the other locked windows remain unchanged.

`computation_config.json` records input checks, input/source hashes, exact parameters, random streams and actual server dependency versions. `integration_diagnostics.json` records quadrature and independent survival-integral checks. The combined `cox_batches.csv` and `wr_batches.csv` contain a case column and support numerical-error checks. Source and runtime provenance is retained in the configuration; it does not change the manuscript's R 4.4.x notation. `SHA256SUMS.txt` covers every result file except itself.

This is a new H1-observed computation, with new numerical methods and random streams. It does not claim reproduction of the archived H1-standard numbers. Manuscript tables and existing code bundles are unchanged.

Use the submitted population-estimand workflow with its fixed configuration:

```text
python -B workflows/population_estimands/compute_h1_observed_estimands.py --recompute --output-dir /path/to/new_external_results
```

The entry imports the shared numerical implementation from `src/dgm/`. Choose a new output directory for each calculation.
"""
        (output / "README.md").write_text(readme, encoding="utf-8")
        configuration["status"] = "completed"
        configuration["completed_utc"] = utc_now()
        configuration["result_files_sha256"] = {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file() and p != config_path}
        save_config()
        manifest = "".join(f"{sha256(p)}  {p.name}\n" for p in sorted(output.iterdir()) if p.is_file() and p.name != "SHA256SUMS.txt")
        (output / "SHA256SUMS.txt").write_text(manifest, encoding="utf-8")
        print("Completed. Results saved to " + str(output), flush=True)
    except BaseException as error:
        configuration["status"] = "failed"
        configuration["failed_utc"] = utc_now()
        configuration["failure"] = {"type": type(error).__name__, "message": str(error)}
        save_config()
        raise


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_portable_config():
    path = PACKAGE / "config/estimands/H1_observed_verification.json"
    config = json.loads(path.read_text(encoding="utf-8"))
    require(sha256(Path(__file__)) == config["portable_sha256"][Path(__file__).name],
            "Portable entry SHA-256 mismatch.")
    require(sha256(Path(dgm.__file__)) == config["current_helper_sha256"],
            "Locked numerical helper SHA-256 mismatch.")
    require(sha256(PACKAGE / config["scientific_configuration"]["path"]) == config["scientific_configuration"]["sha256"], "Scientific calculation settings differ.")
    return config


def verify_locked_inputs(input_dir, config):
    for name, expected in config["input_sha256"].items():
        require(sha256(input_dir / name) == expected, f"Locked input SHA-256 mismatch: {name}")
    return verify_inputs(input_dir, dgm.build_specs())


def verify_archive(archive, config):
    expected = config["archive_sha256"]
    require(archive.is_dir(), f"Missing result archive; supply --archive-dir: {archive}")
    for name, digest in expected.items():
        require((archive / name).is_file(), f"Missing saved result: {name}")
        require(sha256(archive / name) == digest, f"Saved result SHA-256 mismatch: {name}")
    saved_config = json.loads((PACKAGE / "config/estimands/H1_observed_calculation.json").read_text(encoding="utf-8"))
    require(saved_config["source_sha256"] == config["source_sha256"],
            "Scientific source identities differ from the verification specification.")
    require({("dgm_parameters.csv" if k == "DGM_all_scenarios.csv" else k): v for k, v in saved_config["input_verification"]["input_sha256"].items()} == config["input_sha256"],
            "Scientific calculation inputs differ from supplied locked inputs.")
    rows = read_csv(archive / "H1_observed_true_estimands.csv")
    cases = config["verified_archive_case_order"]
    require([row["case"] for row in rows] == cases, "Expected four unique cases in locked order.")
    cox = read_csv(archive / "cox_batches.csv")
    wr = read_csv(archive / "wr_batches.csv")
    require(len(cox) == 96 and len(wr) == 400, "Unexpected total batch-row count.")
    diagnostics = json.loads((archive / "integration_diagnostics.json").read_text(encoding="utf-8"))
    require(set(diagnostics) == set(cases), "Integration diagnostic cases differ.")
    specs = {s.case: s for s in dgm.build_specs()}
    fields_checked = 0
    def close(actual, expected_value, label):
        nonlocal fields_checked
        require(math.isclose(float(actual), float(expected_value), rel_tol=1e-12, abs_tol=1e-12),
                f"Saved summary mismatch: {label}")
        fields_checked += 1
    for index, row in enumerate(rows):
        case = row["case"]
        spec = specs[case]
        require(row["scenario"] == "H1-observed", "Unexpected saved scenario.")
        require(row["locked_dgm_version_identifier"] == spec.version and row["generator"] == spec.generator,
                f"Locked specification identity mismatch: {case}")
        seed = 2026102300 + index * 100
        close(row["seed"], seed, case + " seed")
        close(row["tau_rmst_months"], spec.trial.tau_rmst, case + " RMST horizon")
        close(row["tau_rmut_months"], spec.trial.tau_rmut, case + " RMUT horizon")
        close(row["rmut_weight_recurrence_state"], spec.trial.w_r, case + " utility weight")
        for batch_rows, count, stream, size_field, size in (
            (cox, COX_BATCHES, 0, "n_per_arm", COX_N_ARM),
            (wr, WR_BATCHES, 1, "pairs", WR_PAIRS),
        ):
            selected = [r for r in batch_rows if r["case"] == case]
            require([int(r["batch"]) for r in selected] == list(range(1, count + 1)),
                    f"Duplicate/missing batch identifiers: {case}")
            for k, r in enumerate(selected):
                require(r["scenario"] == "H1-observed" and int(r["seed"]) == seed and
                        r["spawn_key"] == f"{stream};{k}" and int(r[size_field]) == size,
                        f"Batch random stream or size mismatch: {case}, {stream}, {k}")
        c = [r for r in cox if r["case"] == case]
        w = [r for r in wr if r["case"] == case]
        logs = [float(r["cox_log_hr"]) for r in c]
        for r in c:
            close(math.log(float(r["cox_hr"])), r["cox_log_hr"], case + " batch log HR")
        log_hr = statistics.mean(logs)
        log_mcse = statistics.stdev(logs) / math.sqrt(COX_BATCHES)
        close(row["cox_log_hr"], log_hr, case + " log HR mean")
        close(row["cox_hr"], math.exp(log_hr), case + " HR projection")
        close(row["cox_log_hr_mcse"], log_mcse, case + " log HR MCSE")
        close(row["cox_hr_mcse"], math.exp(log_hr) * log_mcse, case + " HR MCSE")
        wins = [float(r["conditional_weighted_wins_mean"]) for r in w]
        losses = [float(r["conditional_weighted_losses_mean"]) for r in w]
        ratio = statistics.mean(wins) / statistics.mean(losses)
        se = statistics.stdev([a - ratio * b for a, b in zip(wins, losses)]) / math.sqrt(WR_BATCHES) / statistics.mean(losses)
        close(row["wr"], ratio, case + " WR")
        close(row["wr_mcse"], se, case + " WR MCSE")
        close(row["wr_log"], math.log(ratio), case + " log WR")
        close(row["wr_expected_ipcw_wins_per_pair"], statistics.mean(wins), case + " wins")
        close(row["wr_expected_ipcw_losses_per_pair"], statistics.mean(losses), case + " losses")
        close(row["delta_rmst_months"], float(row["observed_composite_rmst_months"]) - float(row["control_composite_rmst_months"]), case + " RMST contrast")
        close(row["delta_integral_sd_months"], float(row["observed_death_integral_at_rmut_tau_months"]) - float(row["control_death_integral_at_rmut_tau_months"]), case + " death integral contrast")
        close(row["delta_integral_sc_months"], float(row["observed_composite_integral_at_rmut_tau_months"]) - float(row["control_composite_integral_at_rmut_tau_months"]), case + " composite integral contrast")
        weight = spec.trial.w_r
        close(row["delta_rmut_months"], (1 - weight) * float(row["delta_integral_sd_months"]) + weight * float(row["delta_integral_sc_months"]), case + " RMUT contrast")
        require(diagnostics[case]["independent_survival_integral_crosscheck"]["status"] == "passed",
                f"Original independent integration check failed: {case}")
    return {"status": "passed", "saved_cases": 4, "cox_batch_rows": len(cox),
            "wr_batch_rows": len(wr), "saved_files_sha256_verified": len(expected),
            "summary_fields_checked_from_batches_and_integrals": fields_checked,
            "note": "Archive integrity and saved arithmetic were checked. No new numerical integration or Monte Carlo sampling was performed."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--verify-archive", action="store_true", help="Default: verify the existing result archive without NumPy/SciPy.")
    action.add_argument("--recompute", action="store_true", help="Explicitly request the numerical calculation.")
    parser.add_argument("--input-dir", type=Path, default=PACKAGE / "config")
    parser.add_argument("--archive-dir", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--case", choices=[s.case for s in dgm.build_specs()])
    parser.add_argument("--integrals-only", action="store_true", help="With --recompute, calculate only deterministic RMST/RMUT; this is a partial diagnostic, not the full four-estimand table.")
    args = parser.parse_args()
    config = load_portable_config()
    input_dir = args.input_dir.resolve()
    verified = verify_locked_inputs(input_dir, config)
    if not args.recompute:
        require(args.output_dir is None and not args.integrals_only and args.case is None,
                "Output, --case and --integrals-only require --recompute.")
        print(json.dumps({"input_verification": verified, "archive_verification": verify_archive(args.archive_dir.resolve(), config)}, indent=2))
        return
    require(args.output_dir is not None, "Explicit --recompute requires --output-dir.")
    output = args.output_dir.resolve()
    require(not output.exists(), "The output directory must be new.")
    forbidden = [PACKAGE]
    if (PACKAGE.parent / "04_Data_and_Results").is_dir():
        forbidden.append(PACKAGE.parent)
    require(not any(output == p or p in output.parents for p in forbidden),
            "Choose an output directory outside the entire submission/code package.")
    load_numerical_runtime()
    if not args.integrals_only:
        argv = ["--input-dir", str(input_dir), "--output-dir", str(output)]
        if args.case:
            argv += ["--cases", args.case]
        historical_recompute(argv)
    else:
        output.mkdir(parents=True)
        rows, diagnostics = [], {}
        for spec in dgm.build_specs():
            if args.case and spec.case != args.case:
                continue
            row, detail = deterministic_estimands(spec)
            rows.append({"case": spec.case, "scenario": "H1-observed", **row})
            diagnostics[spec.case] = detail
        write_csv(output / "deterministic_RMST_RMUT.csv", rows)
        (output / "integration_diagnostics.json").write_text(json.dumps(diagnostics, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    provenance = {"mode": "deterministic_integrals_only" if args.integrals_only else "full_recompute",
                  "portable_configuration": config, "input_verification": verified,
                  "actual_runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
                  "note": "Saved archives are never overwritten. Numerical errors describe computation, not original-trial or DGM-fit uncertainty."}
    (output / "portable_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    manifest = "".join(f"{sha256(p)}  {p.name}\n" for p in sorted(output.iterdir())
                       if p.is_file() and p.name != "SHA256SUMS.txt")
    (output / "SHA256SUMS.txt").write_text(manifest, encoding="utf-8")
    print("Completed explicit calculation: " + str(output))


if __name__ == "__main__":
    main()
