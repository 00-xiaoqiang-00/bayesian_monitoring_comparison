"""Joint PFS/OS curve calibration for the KEYNOTE-598 baseline DGM.

The numerical model, grids, initializations, bounds and candidate-selection
rule reproduce the project's regularized shared-component specification.
This command fits parameters; it never overwrites the locked configuration.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import time
from dataclasses import asdict, dataclass
from pathlib import Path
import sys
_SUPPORT_CODE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SUPPORT_CODE_ROOT / 'src/common'))
from output_paths import require_external_output

import numpy as np
import scipy
from scipy.optimize import least_squares as _least_squares, minimize
from scipy.stats import norm

REPO_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
INPUT_DIR = REPO_ROOT / "config/dgm_fitting/keynote598"
CURVE_DIR = REPO_ROOT / "inputs/clinical_data/KEYNOTE-598"
RHO = 0.5
MAX_D_SIGMA = 1.5
MIN_R_SIGMA = 0.15
MAX_ERROR_LIMIT = 0.05
MEDIAN_TOLERANCE_MONTHS = 0.25
SURVIVAL_12_TOLERANCE = 0.012
FEASIBILITY_TOLERANCE = 2e-5


def least_squares(*args, **kwargs):
    """Report optimizer progress without changing the optimizer settings."""
    started = time.perf_counter()
    result = _least_squares(*args, **kwargs)
    print(f"  local candidate: nfev={result.nfev}, seconds={time.perf_counter()-started:.2f}", flush=True)
    return result


RISK_TIMES = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33], dtype=float)


RISK = {
    ("os", "treatment"): np.array([284, 245, 223, 204, 180, 146, 100, 62, 37, 10, 0, 0]),
    ("os", "control"): np.array([284, 252, 230, 215, 192, 154, 111, 77, 41, 15, 2, 0]),
    ("pfs", "treatment"): np.array([284, 189, 146, 117, 99, 65, 36, 19, 9, 1, 0, 0]),
    ("pfs", "control"): np.array([284, 198, 157, 124, 104, 75, 43, 20, 11, 4, 1, 0]),
}


EPS = 1e-12


GAUSS_NODES, GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(48)


@dataclass(frozen=True)
class Curve:
    time: np.ndarray
    survival: np.ndarray


def read_curve(endpoint: str, arm: str) -> Curve:
    path = CURVE_DIR / f"KEYNOTE598_{endpoint}_{arm}.csv"
    data = np.loadtxt(path, delimiter=",")
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError(f"Expected two headerless columns: {path}")
    order = np.argsort(data[:, 0], kind="stable")
    time = data[order, 0].astype(float)
    survival = np.minimum.accumulate(np.clip(data[order, 1].astype(float), 0.0, 1.0))
    if time[0] > 0:
        time = np.insert(time, 0, 0.0)
        survival = np.insert(survival, 0, 1.0)
    return Curve(time=time, survival=survival)


def step_value(curve: Curve, t: np.ndarray | float) -> np.ndarray:
    x = np.asarray(t, dtype=float)
    idx = np.searchsorted(curve.time, x, side="right") - 1
    idx = np.clip(idx, 0, len(curve.time) - 1)
    return curve.survival[idx]


def survival_median(time: np.ndarray, survival: np.ndarray) -> float:
    below = np.flatnonzero(survival <= 0.5)
    if len(below) == 0:
        return math.inf
    i = int(below[0])
    if i == 0:
        return float(time[0])
    if survival[i] == survival[i - 1]:
        return float(time[i])
    return float(
        np.interp(
            0.5,
            [survival[i], survival[i - 1]],
            [time[i], time[i - 1]],
        )
    )


def lognormal_survival(t: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    out = np.ones_like(t)
    positive = t > 0
    out[positive] = norm.sf((np.log(t[positive]) - mu) / sigma)
    return np.clip(out, 0.0, 1.0)


def component_pfs_survival(
    t: np.ndarray,
    mu_r: float,
    sigma_r: float,
    mu_d: float,
    sigma_d: float,
    rho: float,
) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    out = np.ones_like(t)
    positive = t > 0
    log_t = np.log(t[positive])
    a = (log_t - mu_r) / sigma_r
    b = (log_t - mu_d) / sigma_d
    u0 = norm.cdf(b)
    half_width = 0.5 * (1.0 - u0)
    midpoint = 0.5 * (1.0 + u0)
    u = midpoint[:, None] + half_width[:, None] * GAUSS_NODES[None, :]
    z_d = norm.ppf(np.clip(u, EPS, 1.0 - EPS))
    conditional = norm.sf((a[:, None] - rho * z_d) / math.sqrt(1.0 - rho**2))
    out[positive] = half_width * (conditional @ GAUSS_WEIGHTS)
    return np.clip(out, 0.0, 1.0)


def curve_fit_grid(curve: Curve, nrisk: np.ndarray, points: int = 520) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    first_zero = int(np.flatnonzero(nrisk == 0)[0])
    end = min(float(curve.time[-1]), float(RISK_TIMES[first_zero]))
    grid = np.linspace(0.03, end - 1e-4, points)
    target = np.interp(grid, curve.time, curve.survival, left=1.0, right=curve.survival[-1])
    risk_fraction = np.interp(grid, RISK_TIMES, nrisk, left=nrisk[0], right=0) / nrisk[0]
    weight = np.sqrt(np.maximum(risk_fraction, 0.005))
    return grid, target, weight


@dataclass(frozen=True)
class SharedMixture:
    pi: float
    mu_r1: float
    sigma_r1: float
    mu_r2: float
    sigma_r2: float
    mu_d1: float
    sigma_d1: float
    mu_d2: float
    sigma_d2: float


def unpack(theta: np.ndarray) -> SharedMixture:
    (
        logit_pi,
        mu_r1,
        log_sigma_r1,
        log_delta_r,
        log_sigma_r2,
        mu_d1,
        log_sigma_d1,
        log_delta_d,
        log_sigma_d2,
    ) = map(float, theta)
    pi = 1.0 / (1.0 + math.exp(-logit_pi))
    return SharedMixture(
        pi=pi,
        mu_r1=mu_r1,
        sigma_r1=math.exp(log_sigma_r1),
        mu_r2=mu_r1 + math.exp(log_delta_r),
        sigma_r2=math.exp(log_sigma_r2),
        mu_d1=mu_d1,
        sigma_d1=math.exp(log_sigma_d1),
        mu_d2=mu_d1 + math.exp(log_delta_d),
        sigma_d2=math.exp(log_sigma_d2),
    )


def pack(
    pi: float,
    mu_r1: float,
    sigma_r1: float,
    mu_r2: float,
    sigma_r2: float,
    mu_d1: float,
    sigma_d1: float,
    mu_d2: float,
    sigma_d2: float,
) -> np.ndarray:
    pi = min(max(float(pi), 0.050001), 0.949999)
    return np.array(
        [
            math.log(pi / (1.0 - pi)),
            mu_r1,
            math.log(sigma_r1),
            math.log(max(mu_r2 - mu_r1, 0.100001)),
            math.log(sigma_r2),
            mu_d1,
            math.log(sigma_d1),
            math.log(max(mu_d2 - mu_d1, 0.100001)),
            math.log(sigma_d2),
        ]
    )


def pfs_survival(t: np.ndarray, fit: SharedMixture) -> np.ndarray:
    s1 = component_pfs_survival(
        t,
        fit.mu_r1,
        fit.sigma_r1,
        fit.mu_d1,
        fit.sigma_d1,
        RHO,
    )
    s2 = component_pfs_survival(
        t,
        fit.mu_r2,
        fit.sigma_r2,
        fit.mu_d2,
        fit.sigma_d2,
        RHO,
    )
    return np.clip(fit.pi * s1 + (1.0 - fit.pi) * s2, 0.0, 1.0)


def os_survival(t: np.ndarray, fit: SharedMixture) -> np.ndarray:
    s1 = lognormal_survival(t, fit.mu_d1, fit.sigma_d1)
    s2 = lognormal_survival(t, fit.mu_d2, fit.sigma_d2)
    return np.clip(fit.pi * s1 + (1.0 - fit.pi) * s2, 0.0, 1.0)


def starts_for_arm(
    arm: str,
    current: dict[str, dict[str, float]],
    oak: dict[str, dict[str, float]],
) -> list[np.ndarray]:
    c = current[arm]
    o = oak[arm]
    averages = [
        0.5 * (c["pi_R1"] + c["pi_D1"]),
        o["pi_R1"],
        c["pi_R1"],
        c["pi_D1"],
    ]
    starts = [
        pack(
            pi,
            o["mu_R1"],
            o["sigma_R1"],
            o["mu_R2"],
            o["sigma_R2"],
            c["mu_D1"],
            c["sigma_D1"],
            c["mu_D2"],
            c["sigma_D2"],
        )
        for pi in averages
    ]
    starts.extend(
        [
            pack(
                0.20,
                math.log(2.0),
                0.6,
                math.log(22.0),
                0.9,
                math.log(1.5),
                0.9,
                math.log(27.0),
                1.1,
            ),
            pack(
                0.35,
                math.log(3.0),
                0.8,
                math.log(25.0),
                1.2,
                math.log(3.0),
                1.2,
                math.log(30.0),
                1.0,
            ),
        ]
    )
    return starts


def coarse_constraints(
    fit: SharedMixture,
    pfs_curve: Curve,
    os_curve: Curve,
    arm: str,
) -> np.ndarray:
    def one(
        model,
        curve: Curve,
        endpoint: str,
    ) -> np.ndarray:
        nrisk = RISK[(endpoint, arm)]
        first_zero = int(np.flatnonzero(nrisk == 0)[0])
        end = min(float(curve.time[-1]), float(RISK_TIMES[first_zero]))
        grid = np.linspace(0.03, end - 1e-4, 500)
        risk = np.interp(
            grid, RISK_TIMES, nrisk, left=nrisk[0], right=0
        )
        grid = grid[risk >= 5]
        target = np.interp(grid, curve.time, curve.survival)
        delta = model(grid) - target
        median = survival_median(curve.time, curve.survival)
        median_values = model(
            np.array(
                [
                    max(0.001, median - MEDIAN_TOLERANCE_MONTHS),
                    median + MEDIAN_TOLERANCE_MONTHS,
                ]
            )
        )
        target_12 = float(step_value(curve, 12.0))
        fitted_12 = float(model(np.array([12.0]))[0])
        return np.concatenate(
            (
                MAX_ERROR_LIMIT - np.abs(delta),
                np.array(
                    [
                        median_values[0] - 0.5,
                        0.5 - median_values[1],
                        SURVIVAL_12_TOLERANCE
                        - abs(fitted_12 - target_12),
                    ]
                ),
            )
        )

    return np.concatenate(
        (
            one(lambda t: pfs_survival(t, fit), pfs_curve, "pfs"),
            one(lambda t: os_survival(t, fit), os_curve, "os"),
        )
    )


def fit_arm(
    arm: str,
    pfs_curve: Curve,
    os_curve: Curve,
    starts: list[np.ndarray],
) -> tuple[SharedMixture, dict[str, object]]:
    pfs_grid, pfs_target, pfs_weight = curve_fit_grid(
        pfs_curve, RISK[("pfs", arm)]
    )
    os_grid, os_target, os_weight = curve_fit_grid(
        os_curve, RISK[("os", arm)]
    )
    pfs_median = survival_median(pfs_curve.time, pfs_curve.survival)
    os_median = survival_median(os_curve.time, os_curve.survival)
    pfs_12 = float(step_value(pfs_curve, 12.0))
    os_12 = float(step_value(os_curve, 12.0))

    def residual(theta: np.ndarray) -> np.ndarray:
        fit = unpack(theta)
        anchors = np.array(
            [
                (float(pfs_survival(np.array([pfs_median]), fit)[0]) - 0.5)
                / 0.008,
                (float(os_survival(np.array([os_median]), fit)[0]) - 0.5)
                / 0.008,
                (float(pfs_survival(np.array([12.0]), fit)[0]) - pfs_12)
                / 0.008,
                (float(os_survival(np.array([12.0]), fit)[0]) - os_12)
                / 0.008,
            ]
        )
        return np.concatenate(
            (
                (pfs_survival(pfs_grid, fit) - pfs_target)
                * pfs_weight
                / 0.02,
                (os_survival(os_grid, fit) - os_target)
                * os_weight
                / 0.02,
                np.repeat(anchors, 8),
            )
        )

    bounds = (
        np.array(
            [
                math.log(0.05 / 0.95),
                math.log(0.15),
                math.log(MIN_R_SIGMA),
                math.log(0.10),
                math.log(MIN_R_SIGMA),
                math.log(0.10),
                math.log(0.12),
                math.log(0.10),
                math.log(0.12),
            ]
        ),
        np.array(
            [
                math.log(0.95 / 0.05),
                math.log(100.0),
                math.log(3.0),
                math.log(300.0),
                math.log(3.0),
                math.log(100.0),
                math.log(MAX_D_SIGMA),
                math.log(300.0),
                math.log(MAX_D_SIGMA),
            ]
        ),
    )
    local_candidates = [
        least_squares(
            residual,
            np.clip(start, bounds[0] + 1e-6, bounds[1] - 1e-6),
            bounds=bounds,
            max_nfev=18000,
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
        )
        for start in starts
    ]
    best = min(local_candidates, key=lambda x: float(np.dot(x.fun, x.fun)))

    def objective(theta: np.ndarray) -> float:
        r = residual(theta)
        return float(np.dot(r, r) / len(r))

    def constraints(theta: np.ndarray) -> np.ndarray:
        return coarse_constraints(unpack(theta), pfs_curve, os_curve, arm)

    constrained = minimize(
        objective,
        best.x,
        method="SLSQP",
        bounds=list(zip(bounds[0], bounds[1])),
        constraints={"type": "ineq", "fun": constraints},
        options={"maxiter": 700, "ftol": 1e-12, "disp": False},
    )
    candidates = [best, constrained]
    feasible = [
        x
        for x in candidates
        if np.min(constraints(x.x)) >= -FEASIBILITY_TOLERANCE
    ]
    selected = min(
        feasible if feasible else candidates,
        key=lambda x: (
            max(0.0, -float(np.min(constraints(x.x)))),
            float(objective(x.x)),
        ),
    )
    info = {
        "least_squares_candidates": len(local_candidates),
        "constrained_success": bool(constrained.success),
        "constrained_message": str(constrained.message),
        "selected_feasible_on_constraint_grid": bool(
            np.min(constraints(selected.x)) >= -FEASIBILITY_TOLERANCE
        ),
        "selected_minimum_constraint": float(np.min(constraints(selected.x))),
        "selected_objective": float(objective(selected.x)),
    }
    return unpack(selected.x), info


def read_initializations():
    """Read pre-existing candidate fits used only to form six initial values."""
    with (INPUT_DIR/"initialization_independent_margin.csv").open(encoding="utf-8-sig",newline="") as f:
        current={r["arm"]:{k:float(r[k]) for k in ("pi_R1","mu_R1","sigma_R1","mu_R2","sigma_R2","pi_D1","mu_D1","sigma_D1","mu_D2","sigma_D2")} for r in csv.DictReader(f) if abs(float(r["rho"])-RHO)<=1e-12}
    with (INPUT_DIR/"initialization_os_locked_candidate.csv").open(encoding="utf-8-sig",newline="") as f:
        oak={r["arm"]:{k:float(r[k]) for k in ("pi_R1","mu_R1","sigma_R1","mu_R2","sigma_R2","mu_D","sigma_D")} for r in csv.DictReader(f)}
    if set(current)!={"treatment","control"} or set(oak)!={"treatment","control"}:
        raise ValueError("Initialization files must contain both trial arms")
    for label,table in (("independent-margin",current),("OS-locked candidate",oak)):
        for arm,row in table.items():
            if not all(math.isfinite(v) for v in row.values()):
                raise ValueError(f"Nonfinite initialization in {label}, {arm}")
            if any(v<=0 for k,v in row.items() if k.startswith("sigma_")):
                raise ValueError(f"Nonpositive initialization scale in {label}, {arm}")
            if any(not 0<v<1 for k,v in row.items() if k.startswith("pi_")):
                raise ValueError(f"Invalid initialization weight in {label}, {arm}")
    return current,oak


def read_reference():
    with (INPUT_DIR/"reference_parameters.csv").open(encoding="utf-8-sig",newline="") as f:
        table={r["arm"]:SharedMixture(float(r["pi_1"]),float(r["mu_R1"]),float(r["sigma_R1"]),float(r["mu_R2"]),float(r["sigma_R2"]),float(r["mu_D1"]),float(r["sigma_D1"]),float(r["mu_D2"]),float(r["sigma_D2"])) for r in csv.DictReader(f)}
    if set(table)!={"treatment","control"}:raise ValueError("Reference parameters must contain both arms")
    for arm,fit in table.items():
        if not all(math.isfinite(v) for v in asdict(fit).values()):raise ValueError(f"Nonfinite reference parameter: {arm}")
        if not 0<fit.pi<1 or min(fit.sigma_r1,fit.sigma_r2,fit.sigma_d1,fit.sigma_d2)<=0:raise ValueError(f"Invalid reference weight/scale: {arm}")
    return table


def validate_curve_file(path):
    values=np.loadtxt(path,delimiter=",")
    if values.ndim!=2 or values.shape[1]!=2 or len(values)<2:
        raise ValueError(f"Expected at least two time/survival rows: {path.name}")
    if not np.isfinite(values).all():raise ValueError(f"Nonfinite curve coordinate: {path.name}")
    if np.any(values[:,0]<0) or np.any(np.diff(values[:,0])<0):raise ValueError(f"Curve times must be nonnegative and ordered: {path.name}")
    if np.any((values[:,1]<0)|(values[:,1]>1)):raise ValueError(f"Survival probabilities must be in [0,1]: {path.name}")


def write_csv(path, rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def sha256(path):
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    global INPUT_DIR, CURVE_DIR
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir",type=Path,default=INPUT_DIR)
    parser.add_argument("--curve-dir",type=Path,default=CURVE_DIR)
    parser.add_argument("--output-dir",type=Path,required=True,help="New output directory outside the code package")
    parser.add_argument("--arm",choices=("both","treatment","control"),default="both")
    parser.add_argument("--require-reference-match",action="store_true",help="Also require fitted parameters to be within --parameter-tolerance of the recorded reference; never influences estimation")
    parser.add_argument("--parameter-tolerance",type=float,default=1e-4)
    args=parser.parse_args()
    INPUT_DIR=args.input_dir.resolve()
    CURVE_DIR=args.curve_dir.resolve()
    output=require_external_output(args.output_dir, _SUPPORT_CODE_ROOT)
    if output.is_relative_to(REPO_ROOT):
        raise ValueError("Output must be outside the code package")
    if output.is_relative_to(INPUT_DIR) or output.is_relative_to(CURVE_DIR):
        raise ValueError("Output must be outside both custom input directories")
    if output.exists():
        raise FileExistsError("Select a new output directory; existing results are never overwritten")
    required=[CURVE_DIR/f"KEYNOTE598_{endpoint}_{arm}.csv" for endpoint in ("pfs","os") for arm in ("treatment","control")]
    required += [INPUT_DIR/"initialization_independent_margin.csv",INPUT_DIR/"initialization_os_locked_candidate.csv",INPUT_DIR/"reference_parameters.csv"]
    missing=[str(p) for p in required if not p.is_file()]
    if missing:raise FileNotFoundError("Missing inputs: "+", ".join(missing))
    if not math.isfinite(args.parameter_tolerance) or args.parameter_tolerance<=0:raise ValueError("Parameter tolerance must be finite and positive")
    for path in required[:4]:validate_curve_file(path)
    curves={(endpoint,arm):read_curve(endpoint,arm) for endpoint in ("pfs","os") for arm in ("treatment","control")}
    current,oak=read_initializations();reference=read_reference()
    arms=("treatment","control") if args.arm=="both" else (args.arm,)
    output.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter();parameter_rows=[];optimizer_rows=[];comparison_rows=[];curve_rows=[]
    for arm in arms:
        print(f"Fitting {arm}: shared weight and both endpoint component parameters",flush=True)
        fit,info=fit_arm(arm,curves[("pfs",arm)],curves[("os",arm)],starts_for_arm(arm,current,oak))
        optimizer_rows.append({"arm":arm,**info})
        parameter_rows.append({"arm":arm,"rho":RHO,**asdict(fit)})
        for name,value in asdict(fit).items():
            expected=getattr(reference[arm],name)
            comparison_rows.append({"arm":arm,"parameter":name,"refitted":value,"reference":expected,"absolute_difference":abs(value-expected),"within_parameter_tolerance":abs(value-expected)<=args.parameter_tolerance})
        grid=np.linspace(0.0,33.0,1000)
        for endpoint,fun in (("pfs",pfs_survival),("os",os_survival)):
            fitted=fun(grid,fit);ref=fun(grid,reference[arm])
            for t,value,expected in zip(grid,fitted,ref):
                curve_rows.append({"arm":arm,"endpoint":endpoint,"time_months":t,"refitted_survival":value,"reference_survival":expected,"absolute_difference":abs(value-expected)})
        print(f"Finished {arm}: pi={fit.pi:.12f}, feasible={info['selected_feasible_on_constraint_grid']}",flush=True)
    write_csv(output/"fitted_parameters.csv",parameter_rows)
    write_csv(output/"optimizer_status.csv",optimizer_rows)
    write_csv(output/"reference_parameter_comparison.csv",comparison_rows)
    write_csv(output/"reference_curve_comparison.csv",curve_rows)
    feasible=all(r["selected_feasible_on_constraint_grid"] for r in optimizer_rows)
    reference_match=all(r["within_parameter_tolerance"] for r in comparison_rows)
    report={"case":"KEYNOTE-598","model":"two-component shared-label bivariate lognormal","calibration":"joint PFS and OS marginal-curve matching","rho":RHO,"max_death_log_sd":MAX_D_SIGMA,"min_progression_log_sd":MIN_R_SIGMA,"all_selected_feasible":feasible,"reference_parameter_tolerance":args.parameter_tolerance,"all_parameters_match_reference":reference_match,"max_reference_parameter_difference":max(r["absolute_difference"] for r in comparison_rows),"max_reference_curve_difference":max(r["absolute_difference"] for r in curve_rows),"elapsed_seconds":time.perf_counter()-started,"runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},"input_sha256":{p.name:sha256(p) for p in required},"optimizer_status":optimizer_rows}
    (output/"fit_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in ("input_sha256","optimizer_status")},indent=2),flush=True)
    return 0 if feasible and (not args.require_reference_match or reference_match) else 2


if __name__=="__main__":
    raise SystemExit(main())
