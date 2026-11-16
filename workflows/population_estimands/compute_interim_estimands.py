"""Reproduce the scheduled DGM RMST/RMUT contrasts reported in Table S2E.

The calculation retains the fixed 1,500-interval trapezoidal
integration and bivariate-normal survival evaluation. It generates no patients,
fits no models, and does not rerun the trial-monitoring simulations.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

PACKAGE = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
CASES = ("survival_colon", "OAK", "IMpower131", "KEYNOTE-598")
SCENARIOS = ("H1-standard", "H1-observed")
SCHEDULE_NAME = "analysis_schedule.csv"
FIELDS = ("case", "analysis", "calendar_m", "tau_rmst_months", "tau_rmut_months",
          "scenario", "dRMST", "dRMUT")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_inputs(package: Path, config: dict) -> tuple[list, dict, dict, dict]:
    folder = package / "config"
    hashes = {}
    for name, expected in config["input_sha256"].items():
        path = folder / name
        require(path.is_file(), f"Missing locked input: {path}")
        hashes[name] = sha256(path)
        require(hashes[name] == expected, f"Locked input SHA-256 mismatch: {name}")
    parameters = read_csv(folder / "dgm_parameters.csv")
    schedule = read_csv(folder / SCHEDULE_NAME)
    trial_rows = read_csv(folder / "trial_process.csv")
    trials = {row["case"]: row for row in trial_rows}
    require(len(trials) == len(trial_rows) == 4 and set(trials) == set(CASES),
            "Expected one trial-process row for each of the four locked cases.")
    analyses = {}
    for case in CASES:
        selected = [row for row in schedule if row["case"] == case and
                    (row["analysis"].startswith("主分析") or row["analysis"] == "最终分析")]
        require(len(selected) == 4 and len({r["analysis"] for r in selected}) == 4,
                f"Expected three primary interim looks and a final look: {case}")
        interim = [r for r in selected if r["analysis"].startswith("主分析")]
        final = [r for r in selected if r["analysis"] == "最终分析"]
        require(len(interim) == 3 and len(final) == 1, f"Invalid analysis schedule: {case}")
        require([float(r["calendar_months"]) for r in selected] ==
                sorted(float(r["calendar_months"]) for r in selected),
                f"Primary schedule is not chronological: {case}")
        require(selected[-1]["analysis"] == "最终分析", f"Final look is not last: {case}")
        for row in selected:
            require(all(math.isfinite(float(row[k])) and float(row[k]) > 0
                        for k in ("RMST_tau_months", "RMUT_tau_months")),
                    f"Invalid restriction window: {case}, {row['analysis']}")
        weight = float(trials[case]["RMUT_w_R"])
        require(math.isfinite(weight) and 0 <= weight <= 1, f"Invalid RMUT weight: {case}")
        analyses[case] = selected
        for scenario in SCENARIOS:
            for arm in ("treatment", "control"):
                components(parameters, case, scenario, arm)
    return parameters, analyses, trials, hashes


def components(rows: list, case: str, scenario: str, arm: str) -> list:
    selected = [r for r in rows if r["case"] == case and
                r["scenario"] == scenario and r["arm"] == arm]
    require(bool(selected), f"Missing DGM components: {case}, {scenario}, {arm}")
    require(len({r["component"] for r in selected}) == len(selected), "Duplicate DGM components.")
    require(math.isclose(sum(float(r["weight"]) for r in selected), 1.0, abs_tol=1e-9),
            f"Mixture weights do not sum to one: {case}, {scenario}, {arm}")
    for row in selected:
        for key in ("weight", "mu_log_R", "sigma_log_R", "mu_log_D", "sigma_log_D", "rho_log_R_D"):
            require(math.isfinite(float(row[key])), f"Non-finite DGM field: {key}")
        require(float(row["weight"]) >= 0 and float(row["sigma_log_R"]) > 0 and
                float(row["sigma_log_D"]) > 0 and -1 < float(row["rho_log_R_D"]) < 1,
                f"Invalid bivariate-lognormal component: {case}, {scenario}, {arm}")
    # Preserve CSV component order, as in the calculation.
    return selected


def expected_keys(analyses: dict) -> list[tuple[str, str, str]]:
    return [(case, look["analysis"], scenario) for case in CASES
            for look in analyses[case] for scenario in SCENARIOS]


def validate_result(rows: list, analyses: dict) -> dict:
    require(len(rows) == 32, "Expected exactly 32 scheduled contrast rows.")
    actual_keys = [(r["case"], r["analysis"], r["scenario"]) for r in rows]
    require(actual_keys == expected_keys(analyses), "Result keys/order differ from the locked schedule.")
    lookup = {(case, r["analysis"]): r for case in CASES for r in analyses[case]}
    for row in rows:
        look = lookup[row["case"], row["analysis"]]
        for result_field, input_field in (("calendar_m", "calendar_months"),
                                          ("tau_rmst_months", "RMST_tau_months"),
                                          ("tau_rmut_months", "RMUT_tau_months")):
            require(math.isclose(float(row[result_field]), float(look[input_field]),
                                 rel_tol=0.0, abs_tol=1e-12),
                    f"Result/schedule mismatch: {row['case']}, {row['analysis']}, {result_field}")
        require(all(math.isfinite(float(row[k])) for k in ("dRMST", "dRMUT")),
                "Non-finite population contrast.")
    keynote = next(r for r in rows if r["case"] == "KEYNOTE-598" and
                   r["analysis"] == "最终分析")
    require(float(keynote["tau_rmst_months"]) == 12 and float(keynote["tau_rmut_months"]) == 18,
            "KEYNOTE-598 final horizons must be RMST 12 and RMUT 18 months.")
    return {"rows": 32, "cases": 4, "looks_per_case": 4, "scenarios": list(SCENARIOS),
            "KEYNOTE_598_final_horizons_months": {"RMST": 12, "RMUT": 18}}


def curve_survival(rows, time, np, norm, multivariate_normal):
    # This is the fixed formula and operation order.
    time = np.asarray(time, dtype=float)
    log_time = np.log(np.maximum(time, 1e-9))
    composite = np.zeros_like(time)
    death = np.zeros_like(time)
    for row in rows:
        a = (log_time - float(row["mu_log_R"])) / float(row["sigma_log_R"])
        b = (log_time - float(row["mu_log_D"])) / float(row["sigma_log_D"])
        weight = float(row["weight"])
        rho = float(row["rho_log_R_D"])
        death += weight * (1 - norm.cdf(b))
        points = np.column_stack([-a, -b])
        composite += weight * multivariate_normal(mean=[0, 0], cov=[[1, rho], [rho, 1]]).cdf(points)
    require(bool(np.all(np.isfinite(composite))) and bool(np.all(np.isfinite(death))),
            "Non-finite DGM survival probability.")
    return composite, death


def recompute(parameters, analyses, trials, intervals):
    import numpy as np
    import scipy
    from scipy.integrate import trapezoid
    from scipy.stats import multivariate_normal, norm

    results = []
    for case in CASES:
        weight = float(trials[case]["RMUT_w_R"])
        for look in analyses[case]:
            for scenario in SCENARIOS:
                treatment = components(parameters, case, scenario, "treatment")
                control = components(parameters, case, scenario, "control")
                contrasts = {}
                for field, tau_field in (("dRMST", "RMST_tau_months"), ("dRMUT", "RMUT_tau_months")):
                    horizon = float(look[tau_field])
                    grid = np.linspace(0, horizon, intervals + 1)
                    grid[0] = 1e-9
                    sct, sdt = curve_survival(treatment, grid, np, norm, multivariate_normal)
                    scc, sdc = curve_survival(control, grid, np, norm, multivariate_normal)
                    integrand = sct - scc if field == "dRMST" else (
                        (1 - weight) * (sdt - sdc) + weight * (sct - scc))
                    contrasts[field] = float(trapezoid(integrand, grid))
                results.append({"case": case, "analysis": look["analysis"],
                                "calendar_m": look["calendar_months"],
                                "tau_rmst_months": look["RMST_tau_months"],
                                "tau_rmut_months": look["RMUT_tau_months"],
                                "scenario": scenario, **contrasts})
        print("Integrated " + case, flush=True)
    return results, {"numpy": np.__version__, "scipy": scipy.__version__}


def compare_results(actual, reference, tolerance):
    details = []
    for a, r in zip(actual, reference):
        require(all(a[k] == r[k] for k in ("case", "analysis", "scenario")),
                "Computed/reference key mismatch.")
        item = {k: a[k] for k in ("case", "analysis", "scenario")}
        for field in ("dRMST", "dRMUT"):
            item[field + "_absolute_difference_months"] = abs(float(a[field]) - float(r[field]))
            item[field + "_matches_4_decimals"] = round(float(a[field]), 4) == round(float(r[field]), 4)
        details.append(item)
    maxima = {field: max(r[field + "_absolute_difference_months"] for r in details)
              for field in ("dRMST", "dRMUT")}
    passed = all(v <= tolerance for v in maxima.values())
    return {"status": "passed" if passed else "failed", "rows": len(details),
            "absolute_tolerance_months": tolerance, "maximum_absolute_difference_months": maxima,
            "all_values_match_at_4_decimals": all(r[k] for r in details for k in r if k.endswith("_matches_4_decimals")),
            "row_comparisons": details}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--verify-archive", action="store_true", help="Verify the saved table and locked inputs; default mode.")
    mode.add_argument("--recompute", action="store_true", help="Recompute 32 rows using the original 1,500-interval method.")
    parser.add_argument("--package-root", type=Path, default=PACKAGE)
    parser.add_argument("--config", type=Path, help="Optional path to the fixed computation configuration.")
    parser.add_argument("--reference-csv", type=Path, help="Saved 32-row table; auto-detected in the full submission package.")
    parser.add_argument("--output-dir", type=Path, help="New external output directory, required for recomputation.")
    args = parser.parse_args()
    package = args.package_root.resolve()
    config_path = (args.config or package / "config/estimands/interim_computation_config.json").resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    require(sha256(Path(__file__)) == config["portable_script"]["sha256"], "Portable script identity mismatch.")
    parameters, analyses, trials, input_hashes = load_inputs(package, config)
    default_reference = package.parent / config["saved_result"]["path_from_submission_root"]
    reference = args.reference_csv.resolve() if args.reference_csv else (default_reference if default_reference.is_file() else None)
    reference_rows = None
    if reference is not None:
        require(sha256(reference) == config["saved_result"]["sha256"], "Saved-result SHA-256 mismatch.")
        reference_rows = read_csv(reference)
        require(tuple(reference_rows[0]) == FIELDS, "Saved table columns differ from the released two-horizon schema.")
        validate_result(reference_rows, analyses)
    if not args.recompute:
        if reference_rows is None:
            parser.error("The saved CSV is outside this standalone code directory; provide --reference-csv, or use --recompute.")
        print(json.dumps({"status": "passed", "mode": "verify-archive", **validate_result(reference_rows, analyses),
                          "input_sha256": input_hashes, "saved_result_sha256": sha256(reference),
                          "numerical_integration_executed": False}, ensure_ascii=False))
        return
    if args.output_dir is None:
        parser.error("--recompute requires --output-dir.")
    output = args.output_dir.resolve()
    full_submission = (package.parent / "04_Data_and_Results").is_dir()
    protected = package.parent if full_submission else package
    if output.exists() or output == protected or protected in output.parents:
        parser.error("Choose a new output directory outside the submission/code package.")
    require(config["numerical_method"]["trapezoidal_intervals"] == 1500, "Unexpected historical grid size.")
    results, runtime = recompute(parameters, analyses, trials, 1500)
    validate_result(results, analyses)
    comparison = (compare_results(results, reference_rows, config["comparison_absolute_tolerance_months"])
                  if reference_rows is not None else {"status": "not_requested", "reason": "No reference CSV provided/found."})
    output.mkdir(parents=True)
    result_path = output / "interim_true_estimands.csv"
    with result_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, quoting=csv.QUOTE_ALL)
        writer.writeheader(); writer.writerows(results)
    metadata = {"numerical_method": config["numerical_method"], "runtime": runtime,
                "input_sha256": input_hashes, "script_sha256": sha256(Path(__file__)),
                "fixed_configuration_sha256": sha256(config_path),
                "result_sha256": sha256(result_path), "reference_comparison": comparison,
                "RMUT_weights": {case: float(trials[case]["RMUT_w_R"]) for case in CASES},
                "patients_generated": False, "model_fitting_executed": False,
                "archive_overwritten": False}
    (output / "computation_config.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": comparison["status"], "rows": 32, "output": str(output),
                      "maximum_absolute_difference_months": comparison.get("maximum_absolute_difference_months"),
                      "runtime": runtime}, ensure_ascii=False))
    if comparison["status"] == "failed":
        raise SystemExit("Recomputation did not meet the archived numerical tolerance; inspect the external comparison report.")


if __name__ == "__main__":
    main()
