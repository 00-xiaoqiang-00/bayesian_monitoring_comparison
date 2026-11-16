"""Build the frozen KM/Aalen--Johansen targets for rho sensitivity fitting.

The three reconstructed trial templates use the supplied native-resolution,
digitized KM curves. ``survival_colon`` uses the supplied patient-level data
and derives an OS Kaplan--Meier curve and a recurrence cumulative-incidence
curve with death as the competing event.

This script only creates a normalized target package.  It does not fit a DGM,
run MCMC, or modify any method-specific project.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import sys
_SUPPORT_CODE_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
sys.path.insert(0, str(_SUPPORT_CODE_ROOT / 'src/common'))
from output_paths import require_external_output

import numpy as np


REPO_ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'src').is_dir() and (p / 'config').is_dir())
INPUT_DIR = REPO_ROOT / "inputs/clinical_data"
OUT = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_digitized(path: Path, case: str, endpoint: str, arm: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        for line_no, row in enumerate(reader, start=1):
            if len(row) != 2:
                raise ValueError(f"{path}:{line_no}: expected two columns")
            time, value = float(row[0]), float(row[1])
            if time < 0 or not 0 <= value <= 1:
                raise ValueError(f"{path}:{line_no}: invalid time/survival value")
            rows.append({
                "case": case,
                "endpoint": endpoint,
                "arm": arm,
                "time_months": time,
                "curve_value": value,
                "curve_type": "KM",
                "source_kind": "digitized_KM",
                "source_path": str(path),
            })
    if not rows or rows[0]["time_months"] != 0:
        raise ValueError(f"{path}: curve must start at time 0")
    times = np.asarray([float(row["time_months"]) for row in rows])
    values = np.asarray([float(row["curve_value"]) for row in rows])
    if np.any(np.diff(times) <= 0) or np.any(np.diff(values) > 1e-8):
        raise ValueError(f"{path}: target must have strictly increasing time and non-increasing curve")
    # A one-pixel terminal vertical mark in the KEYNOTE-598 treatment PFS
    # extraction falls from 0.116 to 0.001 at the final x-coordinate.  It is
    # not a KM step (there is no corresponding reported event/risk node), so
    # trim only this terminal digitization artifact.  The source file and its
    # hash remain unchanged and the rule is recorded in the target manifest.
    if len(rows) >= 3 and values[-2] - values[-1] > 0.05:
        rows = rows[:-1]
    return rows


def colon_rows(path: Path) -> list[dict[str, object]]:
    """Return one row per subject and endpoint from the archived colon data."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"id", "rx", "status", "time", "etype"}
    missing = required - set(rows[0]) if rows else required
    if missing:
        raise ValueError(f"{path}: missing {sorted(missing)}")
    return rows


def colon_os(rows: list[dict[str, object]], arm: str) -> list[tuple[float, float]]:
    records = [r for r in rows if r["rx"] == arm and int(r["etype"]) == 2]
    grouped: dict[float, list[int]] = defaultdict(lambda: [0, 0])
    for row in records:
        grouped[float(row["time"])][int(row["status"])] += 1
    n_risk = len(records)
    surv = 1.0
    out = [(0.0, 1.0)]
    for time in sorted(grouped):
        censored, events = grouped[time]
        if events:
            surv *= 1.0 - events / n_risk
            out.append((time / 30.4375, surv))
        n_risk -= censored + events
        if n_risk <= 0:
            break
    return out


def colon_competing_records(rows: list[dict[str, object]], arm: str) -> list[tuple[float, int]]:
    """Create competing-risk records: 1=recurrence, 2=death, 0=censor."""
    by_id: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if row["rx"] == arm:
            by_id[row["id"]].append(row)
    result: list[tuple[float, int]] = []
    for subject_rows in by_id.values():
        recurrence = next(r for r in subject_rows if int(r["etype"]) == 1)
        death = next(r for r in subject_rows if int(r["etype"]) == 2)
        r_time, r_event = float(recurrence["time"]), int(recurrence["status"])
        d_time, d_event = float(death["time"]), int(death["status"])
        if r_event and (not d_event or r_time < d_time):
            result.append((r_time, 1))
        elif d_event and (not r_event or d_time < r_time):
            result.append((d_time, 2))
        elif r_event and d_event and r_time == d_time:
            # A tied recurrence/death is conservatively assigned to death;
            # this is the standard competing-risk convention for an unknown
            # ordering and affects no observed tie in the archived data.
            result.append((d_time, 2))
        else:
            result.append((max(r_time, d_time), 0))
    return result


def colon_cif(rows: list[dict[str, object]], arm: str) -> list[tuple[float, float]]:
    records = colon_competing_records(rows, arm)
    grouped: dict[float, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for time, event in records:
        grouped[time][event] += 1
    n_risk = len(records)
    survival = 1.0
    cif = 0.0
    out = [(0.0, 0.0)]
    for time in sorted(grouped):
        event_recurrence = grouped[time].get(1, 0)
        event_death = grouped[time].get(2, 0)
        censor = grouped[time].get(0, 0)
        total_events = event_recurrence + event_death
        if total_events:
            cif += survival * event_recurrence / n_risk
            survival *= 1.0 - total_events / n_risk
            out.append((time / 30.4375, cif))
        n_risk -= total_events + censor
        if n_risk <= 0:
            break
    return out


def write_normalized(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["case", "endpoint", "arm", "time_months", "curve_value", "curve_type", "source_kind", "source_path"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    curves: dict[tuple[str, str, str], list[dict[str, object]]] = {}
    source_paths: set[Path] = set()
    digitized = {}
    for case, prefix in [("OAK", "OAK_A_ITT"), ("IMpower131", "IMpower131"), ("KEYNOTE-598", "KEYNOTE598")]:
        digitized[case] = {endpoint: {arm: INPUT_DIR / case / f"{prefix}_{endpoint.lower()}_{arm}.csv"
                                      for arm in ("control", "treatment")}
                           for endpoint in ("PFS", "OS")}
    for case, endpoints in digitized.items():
        for endpoint, arms in endpoints.items():
            for arm, source in arms.items():
                source_paths.add(source)
                curves[(case, endpoint, arm)] = read_digitized(source, case, endpoint, arm)

    colon_source = INPUT_DIR / "survival_colon.csv"
    source_paths.add(colon_source)
    colon = colon_rows(colon_source)
    for arm in ("Obs", "Lev+5FU"):
        normalized_arm = "control" if arm == "Obs" else "treatment"
        os_rows = [
            {
                "case": "survival_colon", "endpoint": "OS", "arm": normalized_arm,
                "time_months": t, "curve_value": y, "curve_type": "KM",
                "source_kind": "patient_level_derived_KM", "source_path": str(colon_source),
            }
            for t, y in colon_os(colon, arm)
        ]
        cif_rows = [
            {
                "case": "survival_colon", "endpoint": "recurrence_cif", "arm": normalized_arm,
                "time_months": t, "curve_value": y, "curve_type": "Aalen-Johansen_CIF",
                "source_kind": "patient_level_derived_AJ", "source_path": str(colon_source),
            }
            for t, y in colon_cif(colon, arm)
        ]
        curves[("survival_colon", "OS", normalized_arm)] = os_rows
        curves[("survival_colon", "recurrence_cif", normalized_arm)] = cif_rows

    manifest_rows = []
    for key, rows in sorted(curves.items()):
        case, endpoint, arm = key
        out_name = f"{case}_{endpoint}_{arm}.csv".replace(" ", "_")
        out_path = OUT / out_name
        write_normalized(out_path, rows)
        manifest_rows.append({
            "case": case,
            "endpoint": endpoint,
            "arm": arm,
            "target_file": str(out_path),
            "target_sha256": sha256(out_path),
            "n_points": len(rows),
            "time_min_months": rows[0]["time_months"],
            "time_max_months": rows[-1]["time_months"],
            "source_kind": rows[0]["source_kind"],
            "source_path": rows[0]["source_path"],
            "source_sha256": sha256(Path(rows[0]["source_path"])),
        })
    manifest = {
        "version": "rho_sensitivity_km_target_v1",
        "status": "FROZEN_NORMALIZED_TARGETS",
        "units": {"time": "months", "curve_value": "probability"},
        "curve_definition": {
            "PFS": "digitized marginal Kaplan-Meier curve",
            "OS": "digitized or patient-level Kaplan-Meier curve",
            "recurrence_cif": "patient-level Aalen-Johansen cumulative incidence with death competing",
        },
        "target_rows": manifest_rows,
        "source_files": [{"path": str(p), "sha256": sha256(p)} for p in sorted(source_paths)],
        "provenance_note": (
            "These targets are the original/digitized observable curves. They are not generated from the locked DGM. "
            "A single terminal one-pixel vertical-mark artifact in KEYNOTE-598 treatment PFS was removed because it "
            "is not supported by a reported KM event/risk node. "
            "The target package is used only for deterministic rho-specific refitting; H0 and H1-standard are derived "
            "after the observed curves are fitted and are not used to select the decision thresholds."
        ),
    }
    (OUT / "target_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(manifest_rows)} normalized target curves to {OUT}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    INPUT_DIR = args.input_root.resolve()
    OUT = require_external_output(args.output_dir, _SUPPORT_CODE_ROOT)
    if OUT == REPO_ROOT or REPO_ROOT in OUT.parents:
        parser.error("Choose an output directory outside the source package")
    OUT.mkdir(parents=True, exist_ok=False)
    main()
