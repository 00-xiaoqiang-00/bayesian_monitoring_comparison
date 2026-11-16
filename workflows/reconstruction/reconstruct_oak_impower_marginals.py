"""Reconstruct OAK and IMpower131 marginal pseudo-IPD.

The supplied dense KM coordinates, risk tables and reported event counts are
used to allocate events across risk intervals and position event/censoring
times. PFS and OS are reconstructed separately and cannot be paired row by row.
The final tail interval accounts for subjects after the last risk-table node.
"""

from __future__ import annotations

import csv
import argparse
import hashlib
import json
import math
import platform
from dataclasses import dataclass
from pathlib import Path
import sys
_SUPPORT_CODE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SUPPORT_CODE_ROOT / 'src/common'))
from output_paths import require_external_output

import numpy as np
import pandas as pd
import statsmodels
from statsmodels.duration.hazard_regression import PHReg


HERE = Path(__file__).resolve().parent
REPO_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
INPUT_DIR = REPO_ROOT / "inputs/clinical_data"
OUTPUT_DIR = None
EPS = 1e-12


@dataclass(frozen=True)
class Curve:
    time: np.ndarray
    survival: np.ndarray


@dataclass(frozen=True)
class CurveSpec:
    case: str
    endpoint: str
    arm: str
    display_arm: str
    filename: str
    n: int
    risk_times: tuple[float, ...]
    risk_counts: tuple[int, ...]
    total_events: int
    reported_median: float
    landmarks: tuple[float, ...]


OAK_OS_TIMES = tuple(float(x) for x in range(28))
OAK_PFS_TIMES = tuple(float(x) for x in range(25))
IMPOWER_OS_TIMES = tuple(float(x) for x in range(0, 37, 2))
IMPOWER_PFS_TIMES = tuple(float(x) for x in range(0, 29, 2))


SPECS = (
    CurveSpec(
        "OAK", "PFS", "treatment", "Atezolizumab", "OAK_A_ITT_pfs_treatment.csv", 425,
        OAK_PFS_TIMES,
        (425, 395, 243, 190, 181, 139, 128, 119, 111, 99, 92, 80, 75, 64, 59, 53, 51, 49, 45, 29, 22, 12, 9, 7, 2),
        380, 2.8, (6.0, 12.0),
    ),
    CurveSpec(
        "OAK", "PFS", "control", "Docetaxel", "OAK_A_ITT_pfs_control.csv", 425,
        OAK_PFS_TIMES,
        (425, 385, 283, 223, 198, 142, 110, 91, 81, 60, 50, 41, 38, 29, 27, 23, 19, 17, 16, 13, 11, 5, 4, 2, 2),
        375, 4.0, (6.0, 12.0),
    ),
    CurveSpec(
        "OAK", "OS", "treatment", "Atezolizumab", "OAK_A_ITT_os_treatment.csv", 425,
        OAK_OS_TIMES,
        # Figure 2A prints 263 at month 3. The value 363 is verified against
        # FDA CDTL review Figure 1 (PDF page 13; printed page 12).
        # Month 27 prints 1.
        (425, 407, 382, 363, 342, 326, 305, 279, 260, 248, 234, 223, 218, 205, 198, 188, 175, 163, 157, 141, 116, 74, 54, 41, 28, 15, 4, 1),
        271, 13.8, (12.0, 18.0),
    ),
    CurveSpec(
        "OAK", "OS", "control", "Docetaxel", "OAK_A_ITT_os_control.csv", 425,
        tuple(float(x) for x in range(27)),
        # Figure 2A leaves the month-27 entry blank, so the last auditable
        # published count is month 26 (n=3).  The reconstruction adds its own
        # post-table tail interval rather than treating the blank as zero.
        (425, 390, 365, 336, 311, 286, 263, 236, 219, 195, 179, 168, 151, 140, 132, 123, 116, 104, 98, 90, 70, 51, 37, 28, 16, 6, 3),
        298, 9.6, (12.0, 18.0),
    ),
    CurveSpec(
        "IMpower131", "PFS", "treatment", "A+CnP", "IMpower131_pfs_treatment.csv", 343,
        IMPOWER_PFS_TIMES,
        (343, 294, 257, 172, 134, 88, 61, 44, 33, 24, 18, 12, 5, 4, 2),
        270, 6.3, (6.0, 12.0),
    ),
    CurveSpec(
        "IMpower131", "PFS", "control", "CnP", "IMpower131_pfs_control.csv", 340,
        IMPOWER_PFS_TIMES,
        (340, 279, 227, 128, 79, 48, 28, 21, 12, 11, 6, 4, 3, 2, 2),
        289, 5.6, (6.0, 12.0),
    ),
    CurveSpec(
        "IMpower131", "OS", "treatment", "A+CnP", "IMpower131_os_treatment.csv", 343,
        IMPOWER_OS_TIMES,
        (343, 309, 292, 271, 242, 210, 187, 167, 153, 141, 108, 86, 68, 57, 39, 29, 20, 8, 2),
        228, 14.2, (12.0, 24.0),
    ),
    CurveSpec(
        "IMpower131", "OS", "control", "CnP", "IMpower131_os_control.csv", 340,
        IMPOWER_OS_TIMES,
        (340, 316, 293, 266, 229, 209, 184, 157, 138, 116, 85, 67, 49, 36, 22, 16, 11, 6, 2),
        245, 13.5, (12.0, 24.0),
    ),
)


REPORTED_HR = {
    ("OAK", "PFS"): (0.95, 0.82, 1.10, "stratified"),
    ("OAK", "OS"): (0.73, 0.62, 0.87, "stratified"),
    ("IMpower131", "PFS"): (0.71, 0.60, 0.85, "stratified"),
    ("IMpower131", "OS"): (0.88, 0.73, 1.05, "stratified"),
}


def input_path(spec: CurveSpec) -> Path:
    return INPUT_DIR / spec.case / spec.filename


def read_curve(spec: CurveSpec) -> Curve:
    data = np.loadtxt(input_path(spec), delimiter=",")
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError(f"Expected two headerless columns in {input_path(spec)}")
    order = np.argsort(data[:, 0], kind="stable")
    time = data[order, 0].astype(float)
    survival = np.minimum.accumulate(np.clip(data[order, 1].astype(float), 0.0, 1.0))
    keep = np.concatenate(([True], np.diff(time) > 1e-10))
    time, survival = time[keep], survival[keep]
    if time[0] > 0:
        time = np.insert(time, 0, 0.0)
        survival = np.insert(survival, 0, 1.0)
    return Curve(time=time, survival=survival)


def step_value(curve: Curve, t: np.ndarray | float) -> np.ndarray:
    x = np.asarray(t, dtype=float)
    idx = np.searchsorted(curve.time, x, side="right") - 1
    idx = np.clip(idx, 0, len(curve.time) - 1)
    return curve.survival[idx]


def interpolated_median(curve: Curve) -> float:
    below = np.flatnonzero(curve.survival <= 0.5)
    if not len(below):
        return math.inf
    i = int(below[0])
    if i == 0 or curve.survival[i - 1] == curve.survival[i]:
        return float(curve.time[i])
    return float(np.interp(0.5, [curve.survival[i], curve.survival[i - 1]], [curve.time[i], curve.time[i - 1]]))


def km_from_ipd(time: np.ndarray, status: np.ndarray) -> Curve:
    event_times = np.sort(np.unique(time[status == 1]))
    survival = 1.0
    values: list[float] = []
    for t in event_times:
        at_risk = int(np.sum(time >= t - 1e-10))
        events = int(np.sum((np.abs(time - t) < 1e-9) & (status == 1)))
        if at_risk:
            survival *= 1.0 - events / at_risk
        values.append(survival)
    return Curve(np.concatenate(([0.0], event_times)), np.concatenate(([1.0], np.asarray(values))))


def km_median(curve: Curve) -> float:
    below = np.flatnonzero(curve.survival <= 0.5)
    return math.inf if not len(below) else float(curve.time[int(below[0])])


def allocate_integer_events(raw: np.ndarray, capacity: np.ndarray, total: int) -> np.ndarray:
    if int(capacity.sum()) < total:
        raise ValueError(f"Event total {total} exceeds capacity {int(capacity.sum())}")
    allocation = np.clip(np.rint(raw).astype(int), 0, capacity)
    while int(allocation.sum()) < total:
        candidates = np.flatnonzero(allocation < capacity)
        costs = (allocation[candidates] + 1 - raw[candidates]) ** 2 - (allocation[candidates] - raw[candidates]) ** 2
        allocation[candidates[int(np.argmin(costs))]] += 1
    while int(allocation.sum()) > total:
        candidates = np.flatnonzero(allocation > 0)
        costs = (allocation[candidates] - 1 - raw[candidates]) ** 2 - (allocation[candidates] - raw[candidates]) ** 2
        allocation[candidates[int(np.argmin(costs))]] -= 1
    return allocation


def risk_grid_with_tail(spec: CurveSpec, curve: Curve) -> tuple[np.ndarray, np.ndarray, int]:
    times = np.asarray(spec.risk_times, dtype=float)
    counts = np.asarray(spec.risk_counts, dtype=int)
    published_count = len(times)
    if counts[-1] > 0:
        usual_gap = float(np.median(np.diff(times))) if len(times) > 1 else 1.0
        tail_end = max(float(curve.time[-1]) + max(0.10, 0.05 * usual_gap), float(times[-1]) + usual_gap)
        times = np.append(times, tail_end)
        counts = np.append(counts, 0)
    return times, counts, published_count


def interval_event_times(curve: Curve, lo: float, hi: float, count: int) -> np.ndarray:
    if count == 0:
        return np.empty(0)
    mask = (curve.time >= lo) & (curve.time < hi)
    grid_t = np.concatenate(([lo], curve.time[mask], [hi - 1e-6]))
    grid_s = np.asarray(step_value(curve, grid_t), dtype=float)
    start_s = max(float(grid_s[0]), EPS)
    hazard = np.maximum.accumulate(-np.log(np.clip(grid_s / start_s, EPS, 1.0)))
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


def build_ipd(
    curve: Curve,
    risk_times: np.ndarray,
    risk_counts: np.ndarray,
    events_by_interval: np.ndarray,
    align_to_curve: bool,
) -> tuple[np.ndarray, np.ndarray]:
    all_times: list[float] = []
    all_status: list[int] = []
    current_survival = 1.0
    for i, count_raw in enumerate(events_by_interval):
        count = int(count_raw)
        lo, hi = float(risk_times[i]), float(risk_times[i + 1])
        withdrawals = int(risk_counts[i] - risk_counts[i + 1])
        n_censor = withdrawals - count
        censor_times = (
            lo + (np.arange(n_censor) + 0.37) / n_censor * (hi - lo)
            if n_censor else np.empty(0)
        )
        censor_times = np.clip(censor_times, lo + 2e-5, hi - 1e-5)
        initial = interval_event_times(curve, lo, hi, count)
        if not align_to_curve:
            event_times = initial
        else:
            event_times_list: list[float] = []
            previous = lo
            dense = np.unique(np.concatenate((
                np.linspace(lo + 1e-5, hi - 2e-5, 1000),
                curve.time[(curve.time > lo) & (curve.time < hi)],
            )))
            for j in range(count):
                candidates = dense[dense > previous + 1e-6]
                if not len(candidates):
                    candidates = np.asarray([min(hi - 2e-5, previous + 1e-5)])
                censor_before = np.searchsorted(censor_times, candidates, side="left")
                at_risk = np.maximum(int(risk_counts[i]) - j - censor_before, 1)
                next_survival = current_survival * (1.0 - 1.0 / at_risk)
                target = np.asarray(step_value(curve, candidates), dtype=float)
                preferred = initial[j]
                objective = np.abs(next_survival - target) + 1e-5 * np.abs(candidates - preferred) / max(hi - lo, 1e-6)
                chosen_idx = int(np.argmin(objective))
                chosen = float(candidates[chosen_idx])
                current_survival = float(next_survival[chosen_idx])
                event_times_list.append(chosen)
                previous = chosen
            event_times = np.asarray(event_times_list)
        all_times.extend(event_times.tolist())
        all_status.extend([1] * len(event_times))
        all_times.extend(censor_times.tolist())
        all_status.extend([0] * len(censor_times))
    return np.asarray(all_times), np.asarray(all_status, dtype=int)


def reconstruction_score(curve: Curve, time: np.ndarray, status: np.ndarray, end: float) -> float:
    km = km_from_ipd(time, status)
    grid = np.linspace(0.0, min(float(curve.time[-1]), end), 800)
    delta = np.asarray(step_value(km, grid)) - np.asarray(step_value(curve, grid))
    return float(np.sqrt(np.mean(delta**2)) + 0.10 * np.max(np.abs(delta)))


def align_median_within_interval(
    time: np.ndarray,
    status: np.ndarray,
    target: float,
    risk_times: np.ndarray,
) -> np.ndarray:
    km = km_from_ipd(time, status)
    below = np.flatnonzero(km.survival <= 0.5)
    if not len(below) or not math.isfinite(target):
        return time
    k = int(below[0])
    old = float(km.time[k])
    previous = float(km.time[k - 1]) if k > 0 else 0.0
    following = float(km.time[k + 1]) if k + 1 < len(km.time) else math.inf
    interval = int(np.searchsorted(risk_times, old, side="right") - 1)
    interval = max(0, min(interval, len(risk_times) - 2))
    lower = max(previous + 2e-5, float(risk_times[interval]) + 2e-5)
    upper = min(following - 2e-5, float(risk_times[interval + 1]) - 2e-5)
    if lower >= upper or not (float(risk_times[interval]) <= target < float(risk_times[interval + 1])):
        return time
    adjusted = time.copy()
    adjusted[(status == 1) & (np.abs(time - old) < 1e-8)] = float(np.clip(target, lower, upper))
    return adjusted


def reconstruct(spec: CurveSpec, curve: Curve) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    risk_times, risk_counts, published_count = risk_grid_with_tail(spec, curve)
    withdrawals = (risk_counts[:-1] - risk_counts[1:]).astype(int)
    s_start = np.asarray(step_value(curve, risk_times[:-1]), dtype=float)
    s_end = np.asarray(step_value(curve, risk_times[1:] - 1e-5), dtype=float)
    hazard = -np.log(np.clip(s_end / np.maximum(s_start, EPS), EPS, 1.0))
    avg_risk = 0.5 * (risk_counts[:-1] + risk_counts[1:])
    raw = np.clip(hazard * avg_risk, 0, withdrawals)
    capacity = withdrawals.copy()
    supported_capacity = capacity.copy()
    supported_capacity[hazard < 1e-5] = 0
    if int(supported_capacity.sum()) >= spec.total_events:
        capacity = supported_capacity
    allocation = allocate_integer_events(raw, capacity, spec.total_events)

    time, status = build_ipd(curve, risk_times, risk_counts, allocation, align_to_curve=False)
    best_score = reconstruction_score(curve, time, status, float(risk_times[published_count - 1]))
    # A small deterministic integer refinement improves curve agreement while
    # preserving every risk count and the exact event total.
    for _ in range(35):
        best = None
        for i in np.flatnonzero(allocation > 0):
            for j in np.flatnonzero(allocation < capacity):
                if i == j:
                    continue
                candidate = allocation.copy()
                candidate[i] -= 1
                candidate[j] += 1
                c_time, c_status = build_ipd(curve, risk_times, risk_counts, candidate, align_to_curve=False)
                score = reconstruction_score(curve, c_time, c_status, float(risk_times[published_count - 1]))
                if score < best_score - 1e-8:
                    best_score = score
                    best = (candidate, c_time, c_status)
        if best is None:
            break
        allocation, time, status = best

    time, status = build_ipd(curve, risk_times, risk_counts, allocation, align_to_curve=True)
    time = align_median_within_interval(time, status, interpolated_median(curve), risk_times)
    if len(time) != spec.n or int(status.sum()) != spec.total_events:
        raise RuntimeError(f"Invalid dimensions for {spec.case}/{spec.endpoint}/{spec.arm}")
    return time, status, allocation, risk_times, risk_counts, published_count


def write_ipd(spec: CurveSpec, time: np.ndarray, status: np.ndarray) -> Path:
    path = OUTPUT_DIR / f"{spec.case}_{spec.endpoint.lower()}_{spec.arm}_pseudo_ipd_local_20260731.csv"
    order = np.lexsort((-status, time))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["id", "time_months", "status", "case", "endpoint", "arm"])
        for row_id, idx in enumerate(order, 1):
            writer.writerow([row_id, f"{time[idx]:.6f}", int(status[idx]), spec.case, spec.endpoint, spec.arm])
    return path


def curve_error(curve: Curve, km: Curve, max_time: float) -> dict[str, float]:
    grid = np.linspace(0.0, min(float(curve.time[-1]), max_time), 1000)
    delta = np.asarray(step_value(km, grid)) - np.asarray(step_value(curve, grid))
    return {
        "mae": float(np.mean(np.abs(delta))),
        "rmse": float(np.sqrt(np.mean(delta**2))),
        "max_abs_error": float(np.max(np.abs(delta))),
        "evaluation_end_months": float(grid[-1]),
    }


def cox_result(trt: tuple[np.ndarray, np.ndarray], ctrl: tuple[np.ndarray, np.ndarray]) -> dict[str, float]:
    trt_time, trt_status = trt
    ctrl_time, ctrl_status = ctrl
    time = np.concatenate((trt_time, ctrl_time))
    status = np.concatenate((trt_status, ctrl_status))
    arm = np.concatenate((np.ones(len(trt_time)), np.zeros(len(ctrl_time))))
    fit = PHReg(time, arm[:, None], status=status, ties="efron").fit(disp=0)
    beta = float(fit.params[0])
    se = float(fit.bse[0])
    return {
        "reconstructed_unstratified_hr": math.exp(beta),
        "reconstructed_hr_95ci_low": math.exp(beta - 1.959963984540054 * se),
        "reconstructed_hr_95ci_high": math.exp(beta + 1.959963984540054 * se),
        "log_hr_se": se,
    }




def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    curves: dict[tuple[str, str, str], Curve] = {}
    ipd: dict[tuple[str, str, str], tuple[np.ndarray, np.ndarray]] = {}
    risk_rows: list[dict] = []
    endpoint_arm_rows: list[dict] = []
    curve_rows: list[dict] = []
    allocation_record: dict[str, list[int]] = {}

    for spec in SPECS:
        key = (spec.case, spec.endpoint, spec.arm)
        curve = read_curve(spec)
        curves[key] = curve
        time, status, allocation, risk_times, risk_counts, published_count = reconstruct(spec, curve)
        ipd[key] = (time, status)
        write_ipd(spec, time, status)
        allocation_record["_".join(key)] = allocation.tolist()
        reconstructed_risk = np.asarray([np.sum(time >= t - 1e-9) for t in risk_times[:published_count]], dtype=int)
        for t, target, observed in zip(risk_times[:published_count], risk_counts[:published_count], reconstructed_risk):
            note = ""
            if spec.case == "OAK" and spec.endpoint == "OS" and spec.arm == "treatment" and t == 3:
                note = "Published Figure 2A prints 263; 363 is verified against FDA CDTL review Figure 1 (PDF page 13; printed page 12)."
            risk_rows.append({
                "case": spec.case, "endpoint": spec.endpoint, "arm": spec.arm,
                "time_months": t, "target_at_risk": int(target),
                "reconstructed_at_risk": int(observed), "difference": int(observed - target),
                "pass_exact": bool(observed == target), "note": note,
            })
        km = km_from_ipd(time, status)
        errors = curve_error(curve, km, float(risk_times[published_count - 1]))
        curve_rows.append({"case": spec.case, "endpoint": spec.endpoint, "arm": spec.arm, **errors})
        endpoint_arm_rows.append({
            "case": spec.case, "endpoint": spec.endpoint, "arm": spec.arm,
            "display_arm": spec.display_arm, "n": spec.n,
            "target_events": spec.total_events, "reconstructed_events": int(status.sum()),
            "event_difference": int(status.sum() - spec.total_events),
            "reported_median_months": spec.reported_median,
            "digitized_median_months": interpolated_median(curve),
            "pseudo_ipd_km_median_months": km_median(km),
        })
        for landmark in spec.landmarks:
            endpoint_arm_rows[-1][f"digitized_S_{landmark:g}m"] = float(step_value(curve, landmark))
            endpoint_arm_rows[-1][f"pseudo_ipd_S_{landmark:g}m"] = float(step_value(km, landmark))
            endpoint_arm_rows[-1][f"S_{landmark:g}m_difference"] = float(step_value(km, landmark) - step_value(curve, landmark))

    hr_rows: list[dict] = []
    for case in ("OAK", "IMpower131"):
        for endpoint in ("PFS", "OS"):
            result = cox_result(ipd[(case, endpoint, "treatment")], ipd[(case, endpoint, "control")])
            reported, lower, upper, hr_type = REPORTED_HR[(case, endpoint)]
            hr_rows.append({
                "case": case, "endpoint": endpoint,
                "reported_hr": reported, "reported_hr_95ci_low": lower,
                "reported_hr_95ci_high": upper, "reported_hr_type": hr_type,
                **result,
                "hr_point_difference": result["reconstructed_unstratified_hr"] - reported,
                "interpretation_note": "Reconstruction estimate is unstratified; publication estimate is stratified.",
            })

    risk_df = pd.DataFrame(risk_rows)
    arm_df = pd.DataFrame(endpoint_arm_rows)
    curve_df = pd.DataFrame(curve_rows)
    hr_df = pd.DataFrame(hr_rows)
    risk_df.to_csv(OUTPUT_DIR / "risk_table_validation.csv", index=False, encoding="utf-8-sig")
    arm_df.to_csv(OUTPUT_DIR / "event_median_landmark_validation.csv", index=False, encoding="utf-8-sig")
    curve_df.to_csv(OUTPUT_DIR / "curve_reconstruction_error.csv", index=False, encoding="utf-8-sig")
    hr_df.to_csv(OUTPUT_DIR / "hr_validation.csv", index=False, encoding="utf-8-sig")
    (OUTPUT_DIR / "interval_event_allocations.json").write_text(json.dumps(allocation_record, indent=2), encoding="utf-8")

    manifest_rows = []
    for path in sorted([input_path(spec) for spec in SPECS]):
        manifest_rows.append({
            "relative_path": path.relative_to(INPUT_DIR).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": file_sha256(path),
        })
    pd.DataFrame(manifest_rows).to_csv(OUTPUT_DIR / "frozen_input_manifest.csv", index=False, encoding="utf-8-sig")

    summary = {
        "scope": "OAK and IMpower131 marginal pseudo-IPD reconstruction",
        "marginal_warning": "PFS and OS pseudo-IPD are separate marginal reconstructions and cannot be paired row by row.",
        "environment": {
            "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "validation": {
            "all_published_risk_nodes_exact": bool(risk_df["pass_exact"].all()),
            "maximum_absolute_risk_count_error": int(risk_df["difference"].abs().max()),
            "all_total_event_counts_exact": bool((arm_df["event_difference"] == 0).all()),
            "maximum_absolute_event_count_error": int(arm_df["event_difference"].abs().max()),
        },
    }
    (OUTPUT_DIR / "local_reconstruction_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\nHR validation")
    print(hr_df.to_string(index=False))
    print("\nEvent/median/landmark validation")
    print(arm_df.to_string(index=False))
    print("\nCurve errors")
    print(curve_df.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    INPUT_DIR = args.input_root.resolve()
    OUTPUT_DIR = require_external_output(args.output_dir, _SUPPORT_CODE_ROOT)
    if OUTPUT_DIR == REPO_ROOT or REPO_ROOT in OUTPUT_DIR.parents:
        parser.error("Choose an output directory outside the source package.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    main()
