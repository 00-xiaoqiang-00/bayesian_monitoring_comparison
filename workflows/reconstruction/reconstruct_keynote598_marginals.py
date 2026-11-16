"""Reconstruct KEYNOTE-598 marginal pseudo-IPD from fixed curves and risk counts.

PFS and OS are reconstructed separately; their rows are not paired patients.
Mixture-DGM fitting has a separate workflow.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
import sys
_SUPPORT_CODE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SUPPORT_CODE_ROOT / 'src/common'))
from output_paths import require_external_output
import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
INPUT_DIR = REPO_ROOT / "inputs/clinical_data/KEYNOTE-598"
OUTPUT_DIR = None

N_PER_ARM = 284

RISK_TIMES = np.array([0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33], dtype=float)

RISK = {
    ("os", "treatment"): np.array([284, 245, 223, 204, 180, 146, 100, 62, 37, 10, 0, 0]),
    ("os", "control"): np.array([284, 252, 230, 215, 192, 154, 111, 77, 41, 15, 2, 0]),
    ("pfs", "treatment"): np.array([284, 189, 146, 117, 99, 65, 36, 19, 9, 1, 0, 0]),
    ("pfs", "control"): np.array([284, 198, 157, 124, 104, 75, 43, 20, 11, 4, 1, 0]),
}

TOTAL_EVENTS = {
    ("os", "treatment"): 137,
    ("os", "control"): 135,
    ("pfs", "treatment"): 188,
    ("pfs", "control"): 184,
}

REPORTED_MEDIAN = {
    ("os", "treatment"): 21.4,
    ("os", "control"): 21.9,
    ("pfs", "treatment"): 8.2,
    ("pfs", "control"): 8.4,
}

REPORTED_HR = {"os": 1.08, "pfs": 1.06}

EPS = 1e-12

@dataclass(frozen=True)
class Curve:
    time: np.ndarray
    survival: np.ndarray

def read_curve(endpoint: str, arm: str) -> Curve:
    path = INPUT_DIR / f"KEYNOTE598_{endpoint}_{arm}.csv"
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

def km_median(curve: Curve) -> float:
    below = np.flatnonzero(curve.survival <= 0.5)
    return math.inf if len(below) == 0 else float(curve.time[int(below[0])])

def allocate_integer_events(raw: np.ndarray, capacity: np.ndarray, total: int) -> np.ndarray:
    events = np.clip(np.rint(raw).astype(int), 0, capacity)

    while int(events.sum()) < total:
        candidates = np.flatnonzero(events < capacity)
        if len(candidates) == 0:
            raise ValueError("Event total exceeds interval capacity")
        costs = (events[candidates] + 1 - raw[candidates]) ** 2 - (
            events[candidates] - raw[candidates]
        ) ** 2
        events[candidates[int(np.argmin(costs))]] += 1

    while int(events.sum()) > total:
        candidates = np.flatnonzero(events > 0)
        costs = (events[candidates] - 1 - raw[candidates]) ** 2 - (
            events[candidates] - raw[candidates]
        ) ** 2
        events[candidates[int(np.argmin(costs))]] -= 1
    return events

def interval_event_times(curve: Curve, lo: float, hi: float, count: int) -> np.ndarray:
    if count == 0:
        return np.empty(0)
    mask = (curve.time >= lo) & (curve.time < hi)
    grid_t = np.concatenate(([lo], curve.time[mask], [hi - 1e-6]))
    grid_s = np.asarray(step_value(curve, grid_t), dtype=float)
    start_s = max(float(grid_s[0]), EPS)
    hazard = -np.log(np.clip(grid_s / start_s, EPS, 1.0))
    hazard = np.maximum.accumulate(hazard)
    if hazard[-1] <= 1e-10:
        out = lo + (np.arange(count) + 0.5) / count * (hi - lo)
    else:
        unique_h, unique_idx = np.unique(hazard, return_index=True)
        unique_t = grid_t[unique_idx]
        q = (np.arange(count) + 0.5) / count * hazard[-1]
        out = np.interp(q, unique_h, unique_t)
    tiny = min(1e-4, (hi - lo) / (1000 * max(count, 1)))
    out = out + tiny * (np.arange(count) - (count - 1) / 2)
    return np.clip(out, lo + 1e-5, hi - 2e-5)

def build_ipd_from_allocation(
    curve: Curve,
    trisk: np.ndarray,
    nrisk: np.ndarray,
    events_by_interval: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    times: list[float] = []
    status: list[int] = []
    for i, n_events in enumerate(events_by_interval):
        lo, hi = float(trisk[i]), float(trisk[i + 1])
        withdrawals = int(nrisk[i] - nrisk[i + 1])
        n_censor = withdrawals - int(n_events)
        event_times = interval_event_times(curve, lo, hi, int(n_events))
        if n_censor:
            censor_times = lo + (np.arange(n_censor) + 0.37) / n_censor * (hi - lo)
            censor_times = np.clip(censor_times, lo + 2e-5, hi - 1e-5)
        else:
            censor_times = np.empty(0)
        times.extend(event_times.tolist())
        status.extend([1] * len(event_times))
        times.extend(censor_times.tolist())
        status.extend([0] * len(censor_times))
    return np.asarray(times, dtype=float), np.asarray(status, dtype=int)

def build_ipd_aligned_to_curve(
    curve: Curve,
    trisk: np.ndarray,
    nrisk: np.ndarray,
    events_by_interval: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Place fixed event counts at KM-compatible locations within intervals."""
    times: list[float] = []
    status: list[int] = []
    current_survival = 1.0

    for i, n_events_raw in enumerate(events_by_interval):
        n_events = int(n_events_raw)
        lo, hi = float(trisk[i]), float(trisk[i + 1])
        withdrawals = int(nrisk[i] - nrisk[i + 1])
        n_censor = withdrawals - n_events
        if n_censor:
            censor_times = lo + (np.arange(n_censor) + 0.37) / n_censor * (hi - lo)
            censor_times = np.clip(censor_times, lo + 2e-5, hi - 1e-5)
        else:
            censor_times = np.empty(0)

        initial_events = interval_event_times(curve, lo, hi, n_events)
        event_times: list[float] = []
        previous = lo
        dense_candidates = np.unique(
            np.concatenate(
                (
                    np.linspace(lo + 1e-5, hi - 2e-5, 1200),
                    curve.time[(curve.time > lo) & (curve.time < hi)],
                )
            )
        )

        for j in range(n_events):
            candidates = dense_candidates[dense_candidates > previous + 1e-6]
            if len(candidates) == 0:
                candidates = np.array([min(hi - 2e-5, previous + 1e-5)])
            censor_before = np.searchsorted(censor_times, candidates, side="left")
            n_at_risk = np.maximum(int(nrisk[i]) - j - censor_before, 1)
            next_survival = current_survival * (1.0 - 1.0 / n_at_risk)
            target = np.asarray(step_value(curve, candidates), dtype=float)
            preferred = initial_events[j]
            objective = np.abs(next_survival - target) + 1e-5 * np.abs(candidates - preferred) / max(hi - lo, 1e-6)
            chosen_index = int(np.argmin(objective))
            chosen = float(candidates[chosen_index])
            current_survival = float(next_survival[chosen_index])
            event_times.append(chosen)
            previous = chosen

        times.extend(event_times)
        status.extend([1] * len(event_times))
        times.extend(censor_times.tolist())
        status.extend([0] * len(censor_times))

    return np.asarray(times, dtype=float), np.asarray(status, dtype=int)

def km_from_ipd(time: np.ndarray, status: np.ndarray) -> Curve:
    event_times = np.sort(np.unique(time[status == 1]))
    surv = 1.0
    values = []
    for t in event_times:
        n_risk = int(np.sum(time >= t - 1e-10))
        d = int(np.sum((np.abs(time - t) < 1e-9) & (status == 1)))
        if n_risk > 0:
            surv *= 1.0 - d / n_risk
        values.append(surv)
    return Curve(
        time=np.concatenate(([0.0], event_times)),
        survival=np.concatenate(([1.0], np.asarray(values))),
    )

def align_pseudo_ipd_median(time: np.ndarray, status: np.ndarray, target_median: float) -> np.ndarray:
    """Move only the median-crossing event within its risk-table interval."""
    km = km_from_ipd(time, status)
    below = np.flatnonzero(km.survival <= 0.5)
    if len(below) == 0:
        return time
    k = int(below[0])
    old_time = float(km.time[k])
    previous_event = float(km.time[k - 1]) if k > 0 else 0.0
    next_event = float(km.time[k + 1]) if k + 1 < len(km.time) else math.inf
    interval = int(np.searchsorted(RISK_TIMES, old_time, side="right") - 1)
    interval_lo = float(RISK_TIMES[max(interval, 0)])
    interval_hi = float(RISK_TIMES[min(interval + 1, len(RISK_TIMES) - 1)])
    lower = max(previous_event + 2e-5, interval_lo + 2e-5)
    upper = min(next_event - 2e-5, interval_hi - 2e-5)
    if lower >= upper:
        return time
    new_time = float(np.clip(target_median, lower, upper))
    adjusted = time.copy()
    event_mask = (status == 1) & (np.abs(time - old_time) < 1e-8)
    adjusted[event_mask] = new_time
    return adjusted

def reconstruction_score(curve: Curve, time: np.ndarray, status: np.ndarray, nrisk: np.ndarray) -> float:
    km = km_from_ipd(time, status)
    grid = np.linspace(0.0, min(curve.time[-1], RISK_TIMES[np.flatnonzero(nrisk == 0)[0]]), 700)
    target = step_value(curve, grid)
    fitted = step_value(km, grid)
    risk_weight = np.interp(grid, RISK_TIMES, nrisk, left=nrisk[0], right=0) / nrisk[0]
    delta = fitted - target
    return float(np.sqrt(np.average(delta**2, weights=np.maximum(risk_weight, 0.02))) + 0.10 * np.max(np.abs(delta)))

def reconstruct_ipd(endpoint: str, arm: str, curve: Curve) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nrisk_full = RISK[(endpoint, arm)]
    first_zero = int(np.flatnonzero(nrisk_full == 0)[0])
    trisk = RISK_TIMES[: first_zero + 1]
    nrisk = nrisk_full[: first_zero + 1]
    withdrawals = (nrisk[:-1] - nrisk[1:]).astype(int)

    s_start = np.asarray(step_value(curve, trisk[:-1]), dtype=float)
    s_end = np.asarray(step_value(curve, trisk[1:] - 1e-5), dtype=float)
    hazard = -np.log(np.clip(s_end / np.maximum(s_start, EPS), EPS, 1.0))
    average_risk = 0.5 * (nrisk[:-1] + nrisk[1:])
    raw_events = np.clip(hazard * average_risk, 0, withdrawals)
    event_capacity = withdrawals.copy()
    # A flat published tail followed by n.risk=0 indicates censoring, not an
    # artificial final event.  Keep events out of intervals with no visible KM
    # decline and reallocate the published total to supported earlier drops.
    event_capacity[hazard < 1e-4] = 0
    allocation = allocate_integer_events(raw_events, event_capacity, TOTAL_EVENTS[(endpoint, arm)])

    time, status = build_ipd_from_allocation(curve, trisk, nrisk, allocation)
    best_score = reconstruction_score(curve, time, status, nrisk_full)

    # Local integer refinement preserves every risk-table count and the exact
    # total event count while improving agreement with the published curve.
    for _ in range(80):
        best_move = None
        for i in range(len(allocation)):
            if allocation[i] <= 0:
                continue
            for j in range(len(allocation)):
                if i == j or allocation[j] >= event_capacity[j]:
                    continue
                candidate = allocation.copy()
                candidate[i] -= 1
                candidate[j] += 1
                c_time, c_status = build_ipd_from_allocation(curve, trisk, nrisk, candidate)
                score = reconstruction_score(curve, c_time, c_status, nrisk_full)
                if score < best_score - 1e-7:
                    best_score = score
                    best_move = (candidate, c_time, c_status)
        if best_move is None:
            break
        allocation, time, status = best_move

    time, status = build_ipd_aligned_to_curve(curve, trisk, nrisk, allocation)
    time = align_pseudo_ipd_median(
        time,
        status,
        survival_median(curve.time, curve.survival),
    )
    if len(time) != N_PER_ARM or int(status.sum()) != TOTAL_EVENTS[(endpoint, arm)]:
        raise RuntimeError(f"Invalid pseudo-IPD dimensions for {endpoint}/{arm}")
    return time, status, allocation

def risk_counts_from_ipd(time: np.ndarray) -> np.ndarray:
    return np.array([np.sum(time >= t - 1e-9) for t in RISK_TIMES], dtype=int)

def write_ipd(endpoint: str, arm: str, time: np.ndarray, status: np.ndarray) -> Path:
    path = OUTPUT_DIR / f"KEYNOTE598_{endpoint}_{arm}_pseudo_ipd.csv"
    order = np.lexsort((-status, time))
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["id", "time_months", "status", "endpoint", "arm"])
        for row_id, idx in enumerate(order, start=1):
            writer.writerow([row_id, f"{time[idx]:.6f}", int(status[idx]), endpoint.upper(), arm])
    return path

def cox_hr_from_ipd(
    trt_time: np.ndarray,
    trt_status: np.ndarray,
    ctrl_time: np.ndarray,
    ctrl_status: np.ndarray,
) -> float:
    event_times = np.sort(np.unique(np.concatenate((trt_time[trt_status == 1], ctrl_time[ctrl_status == 1]))))

    def score(beta: float) -> float:
        eb = math.exp(beta)
        total = 0.0
        for t in event_times:
            d1 = int(np.sum((np.abs(trt_time - t) < 1e-9) & (trt_status == 1)))
            d0 = int(np.sum((np.abs(ctrl_time - t) < 1e-9) & (ctrl_status == 1)))
            if d1 + d0 == 0:
                continue
            y1 = int(np.sum(trt_time >= t - 1e-10))
            y0 = int(np.sum(ctrl_time >= t - 1e-10))
            total += d1 - (d1 + d0) * y1 * eb / max(y0 + y1 * eb, EPS)
        return total

    lo, hi = -5.0, 5.0
    flo = score(lo)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fm = score(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo = mid
            flo = fm
    return math.exp(0.5 * (lo + hi))

def main() -> None:
    curves = {(endpoint, arm): read_curve(endpoint, arm) for endpoint in ("pfs", "os") for arm in ("treatment", "control")}

    ipd: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    reconstruction_rows = []
    allocations = {}
    for key, curve in curves.items():
        endpoint, arm = key
        time, status, allocation = reconstruct_ipd(endpoint, arm, curve)
        ipd[key] = (time, status)
        allocations[f"{endpoint}_{arm}"] = allocation.tolist()
        write_ipd(endpoint, arm, time, status)
        km = km_from_ipd(time, status)
        nrisk_check = risk_counts_from_ipd(time)
        expected = RISK[key]
        reconstruction_rows.append({
            "endpoint": endpoint.upper(),
            "arm": arm,
            "n": len(time),
            "events": int(status.sum()),
            "reported_events": TOTAL_EVENTS[key],
            "risk_table_max_abs_error": int(np.max(np.abs(nrisk_check - expected))),
            "digitized_median_months": survival_median(curve.time, curve.survival),
            "pseudo_ipd_median_months": km_median(km),
            "reported_median_months": REPORTED_MEDIAN[key],
            "curve_score": reconstruction_score(curve, time, status, expected),
        })

    for endpoint in ("os", "pfs"):
        hr = cox_hr_from_ipd(
            *ipd[(endpoint, "treatment")],
            *ipd[(endpoint, "control")],
        )
        for row in reconstruction_rows:
            if row["endpoint"] == endpoint.upper():
                row["pseudo_ipd_cox_hr"] = hr
                row["reported_hr"] = REPORTED_HR[endpoint]

    with (OUTPUT_DIR / "KEYNOTE598_pseudo_ipd_validation_metrics.csv").open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(reconstruction_rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(reconstruction_rows)

    (OUTPUT_DIR / "KEYNOTE598_pseudo_ipd_interval_event_allocations.json").write_text(
        json.dumps(allocations, indent=2), encoding="ascii"
    )

    print(json.dumps({"case": "KEYNOTE-598", "reconstruction": reconstruction_rows}, indent=2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    INPUT_DIR = args.input_dir.resolve()
    OUTPUT_DIR = require_external_output(args.output_dir, _SUPPORT_CODE_ROOT)
    if OUTPUT_DIR == REPO_ROOT or REPO_ROOT in OUTPUT_DIR.parents:
        parser.error("Choose an output directory outside the source package.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    main()
