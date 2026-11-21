"""Select nine recorded load/speed conditions before any BSE optimization.

The original IK scale and A=+0.1 activation processing are retained. Selection
uses acquisition quality, two angle valleys and CEINMS availability, never a
fitted model's accuracy. Missing torque samples remain explicitly masked.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import find_peaks, lfilter

from audit_subjects import DATA, ROOT, REQUIRED, read_sto
from quality_gate import native_tables, sha, smoothness, gaussian_matlab
from prepare_exp_three_loads import angle_reference, supported
from prepare_subjects import combine, map_predict, structs
from recompute_exp_activation import selected_signal, MODEL_CHANNELS
from review_exp_one_cycle import review_trial
from fit_hill_compliant import input_generator

BASE = ROOT / "results/nine_case_75pct"
VERSION = "excitation-ik-speed-grid-torque-v1"
SUBJECTS = ("240725Exp_A01", "240730Exp_A01")
ADDITIONAL_SMOOTHING_MS = {"240725Exp_A01": 0, "240730Exp_A01": 310}


def candidate_record(subject, entry, choice, model, limits, contact_at_peaks=False):
    """Recover a measured bout with exactly the established signal pipeline."""
    raw_subject = subject.removesuffix("_A01")
    traces, source_label = selected_signal(
        {**entry, "subject": raw_subject, "selected": choice},
        {name: 0.1 for name in MODEL_CHANNELS},
        additional_smoothing_ms=ADDITIONAL_SMOOTHING_MS[subject],
    )
    epoch = traces["timeEpochSeconds"]
    path = DATA / subject / entry["trial"]
    q, angle_audit = angle_reference(path, epoch)
    tables, clocks, identity, exact, issues, repairs = native_tables(path)
    if issues:
        raise ValueError("Native pipeline issues: " + str(issues))
    table = tables["Torque.sto"][1]
    valid = supported(table[:, 0], epoch)
    if not valid.any():
        raise ValueError("No CEINMS torque support in this recorded bout")
    activation = traces["activation"]
    six, front = activation[:, :6], activation[:, 6]
    if (
        not np.isfinite(activation).all()
        or activation.min() < -1e-6
        or activation.max() > 1 + 1e-6
    ):
        raise ValueError("Activation outside [0, 1]")
    back = six @ model["lumpWeights"]
    diagnostic = smoothness(back)
    rough = [key for key, upper in limits.items() if diagnostic[key] > upper]
    if rough:
        raise ValueError("Activation smoothness: " + str(rough))
    torque = np.interp(epoch, table[:, 0], table[:, 1])
    torque[~valid] = np.nan
    load, speed = entry["loadKg"], entry["speedCode"]
    trial_name = Path(entry["trial"]).name
    duration = float(epoch[-1] - epoch[0])
    valleys = find_peaks(-q, prominence=0.12, distance=50)[0]
    correlation = (
        float(np.corrcoef(q[valid], torque[valid])[0, 1]) if valid.sum() > 2 else 0.0
    )
    record = dict(
        id=trial_name + "_nine_case_bout",
        trialPath=entry["trial"],
        loadKg=float(load),
        speedCode=int(speed),
        cadenceBPM=60.0 / duration,
        duration=duration,
        sourceStartSeconds=float(epoch[0] - entry["originEpochSeconds"]),
        sourceEndSeconds=float(epoch[-1] - entry["originEpochSeconds"]),
        q=q,
        torque=torque,
        torqueReferenceValid=valid,
        back=back,
        front=front,
        six=six,
        weights=np.tile(model["lumpWeights"], (len(q), 1)),
        angleTorqueCorrelation=correlation,
        angleUsable=True,
        rightCensored=False,
        split="train",
        fitWeight=1.0,
        boundaryResidualWeight=2.0,
    )
    record["mappedBack"] = map_predict(model, record)
    if ADDITIONAL_SMOOTHING_MS[subject]:
        record["mappedBack"] = gaussian_matlab(record["mappedBack"], 31)
    audit = dict(
        id=record["id"],
        trialPath=record["trialPath"],
        loadKg=float(load),
        speedCode=int(speed),
        sourceStartSeconds=record["sourceStartSeconds"],
        sourceEndSeconds=record["sourceEndSeconds"],
        source=source_label,
        referenceValleys=len(valleys),
        theta=angle_audit,
        torqueReferenceCoverageFraction=float(valid.mean()),
        angleTorqueCorrelation=correlation,
        smoothness=diagnostic,
        torqueIdentity=identity,
        nativeRepairs=repairs,
        maximumSourceGapSeconds=choice["maximumSourceGapSeconds"],
        additionalSubjectWideGaussianWindowMs=ADDITIONAL_SMOOTHING_MS[subject],
        sourceFiles={
            name: sha(path / name) for name in (*REQUIRED, "emg.sto", "emgFilt.sto")
        },
        independentTorqueReferenceAvailable=False,
    )
    if contact_at_peaks:
        contact = np.asarray(choice["peakEpochSeconds"]) - epoch[0]
        if len(contact) != 2 or not 0 < contact[0] < contact[1] < duration:
            raise ValueError("Invalid inferred pick/place times")
        record.update(
            loadContactStartSeconds=float(contact[0]),
            loadContactEndSeconds=float(contact[1]),
        )
        audit.update(
            loadContactStartSeconds=float(contact[0]),
            loadContactEndSeconds=float(contact[1]),
            contactTimingSource="Two recorded EMG peaks; inferred protocol, not measured hand force",
        )
    return record, audit


def select_conditions(candidates):
    """Prefer three distinct speeds per load; never count a repeat as a new case."""
    by_condition = {}
    for record, audit in candidates:
        key = (record["loadKg"], record["speedCode"])
        rank = (
            abs(audit["referenceValleys"] - 2),
            -audit["torqueReferenceCoverageFraction"],
            record["trialPath"],
            record["sourceStartSeconds"],
        )
        if key not in by_condition or rank < by_condition[key][0]:
            by_condition[key] = (rank, record, audit)
    selected = []
    for load in (5, 10, 15):
        options = [v for (mass, _), v in by_condition.items() if mass == load]
        options.sort(
            key=lambda v: (
                abs(v[2]["referenceValleys"] - 2),
                0 if 35 <= v[1]["speedCode"] <= 55 else 1,
                abs(v[1]["speedCode"] - 45),
                -v[2]["torqueReferenceCoverageFraction"],
            )
        )
        selected.extend(options[:3])
    if len(selected) < 9:
        selected_keys = {(v[1]["loadKg"], v[1]["speedCode"]) for v in selected}
        remaining = [v for k, v in by_condition.items() if k not in selected_keys]
        remaining.sort(
            key=lambda v: (
                abs(v[2]["referenceValleys"] - 2),
                -v[2]["torqueReferenceCoverageFraction"],
            )
        )
        selected.extend(remaining[: 9 - len(selected)])
    selected.sort(key=lambda v: (v[1]["loadKg"], v[1]["speedCode"]))
    return [(v[1], v[2]) for v in selected], by_condition


def prepare(subject, output=BASE, audit_only=False, contact_at_peaks=False):
    review = json.loads((ROOT / "docs/exp_one_cycle_review/review.json").read_text())
    old_gate = json.loads(
        (ROOT / "results/exp_subjects_A01/quality_gate.json").read_text()
    )
    limits = old_gate["smoothnessUpperLimits"]
    raw_subject = subject.removesuffix("_A01")
    map_path = ROOT / "results/exp_subjects_A01/map" / (subject + "_map.mat")
    model = loadmat(map_path, simplify_cells=True)["mapModel"]
    cache_path = output / subject / "candidate_review.json"
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    extra_reviews = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    candidates, rejected = [], []
    pattern = re.compile(r"dynamics_stoop_(5|10|15)kg_(25|35|45|55|65)(?:_\d+)?$")
    for path in sorted((DATA / subject / "CalibrationData").iterdir()):
        match = pattern.fullmatch(path.name)
        if not path.is_dir() or not match:
            continue
        trial = str(path.relative_to(DATA / subject))
        entry = next(
            (
                r
                for r in review["conditions"]
                if r["subject"] == raw_subject and r["trial"] == trial
            ),
            None,
        )
        if entry is None:
            if trial not in extra_reviews:
                entry, _ = review_trial(raw_subject, DATA / raw_subject / trial, limits)
                entry.update(loadKg=int(match[1]), speedCode=int(match[2]))
                extra_reviews[trial] = entry
                cache_path.write_text(json.dumps(extra_reviews, indent=2) + "\n")
            entry = extra_reviews[trial]
        seen = set()
        choices = ([entry["selected"]] if entry.get("selected") else []) + entry.get(
            "candidates", []
        )
        for choice in choices:
            key = (
                choice["startEpochSeconds"],
                choice["endEpochSeconds"],
                choice["source"],
            )
            if key in seen or not choice.get("cleanSourceWindow"):
                continue
            seen.add(key)
            try:
                candidates.append(
                    candidate_record(
                        subject, entry, choice, model, limits, contact_at_peaks
                    )
                )
            except (ValueError, IndexError, KeyError) as error:
                rejected.append(
                    dict(trial=trial, startEpochSeconds=key[0], reason=str(error))
                )
    selected, available = select_conditions(candidates)
    records = [r for r, _ in selected]
    gate = dict(
        version=VERSION,
        selectionVersion="nine-recorded-conditions-v1",
        subject=subject,
        requestedConditions=9,
        selectedConditions=len(records),
        uniqueAvailableConditions=len(available),
        selectionUsesFittedPerformance=False,
        allConditionsAreCalibration=True,
        fixedVectorPerSource=True,
        historicalScaleRetained=True,
        missingTorqueIsMasked=True,
        zeroTorqueSupportIsExcluded=True,
        MAPFileSHA256=sha(map_path),
        trials=[a for _, a in selected],
        loadContactAtEMGPeaks=contact_at_peaks,
        loadContactMeasured=False,
        rejectedCandidates=rejected,
    )
    folder = output / subject
    gate_path = folder / "quality_gate.json"
    fitted = (folder / "neural_models.mat").exists() or (
        folder / "checkpoint.mat"
    ).exists()
    if (
        fitted
        and (folder / "prepared.mat").exists()
        and json.loads(gate_path.read_text()) != gate
    ):
        raise RuntimeError(
            "Selection changed after fitting; use a separate result directory"
        )
    gate_path.write_text(json.dumps(gate, indent=2) + "\n")
    for _, a in selected:
        print(
            subject,
            a["loadKg"],
            a["speedCode"],
            a["trialPath"],
            "valleys",
            a["referenceValleys"],
            "torque",
            round(a["torqueReferenceCoverageFraction"], 3),
            flush=True,
        )
    print("AVAILABLE", subject, len(available), "SELECTED", len(records), flush=True)
    if audit_only:
        return
    if len(records) != 9:
        raise RuntimeError(
            "Fewer than nine distinct recorded conditions passed; no duplicates fabricated"
        )

    def generator_mass(r):
        if not contact_at_peaks:
            return np.full(len(r["q"]), r["loadKg"])
        clock = np.arange(len(r["q"])) * 0.01
        raw = r["loadKg"] * (
            (clock >= r["loadContactStartSeconds"])
            & (clock < r["loadContactEndSeconds"])
        )
        decay = np.exp(-0.01 / 0.08)
        return lfilter(1 - decay, [1, -decay], raw)

    generator = input_generator(
        [
            dict(
                a=np.c_[r["back"], r["mappedBack"]],
                front=r["front"],
                m=generator_mass(r),
                loadKg=r["loadKg"],
                q=r["q"],
            )
            for r in records
        ]
    )
    sequences = [combine([r], subject, r["id"]) for r in records]
    # Append continuous three-load demonstrations without resetting between actions.
    by_load = {
        load: [r for r in records if r["loadKg"] == load] for load in (5, 10, 15)
    }
    if all(len(v) == 3 for v in by_load.values()):
        sequences += [
            combine(
                [by_load[load][i] for load in (5, 10, 15)],
                subject,
                "continuous_group_" + str(i + 1),
            )
            for i in range(3)
        ]
    if not fitted or not (folder / "prepared.mat").exists():
        savemat(
            folder / "prepared.mat",
            dict(
                sequences=structs(sequences),
                trainSequences=structs(sequences[:9]),
                cycles=structs(
                    [
                        {k: v for k, v in r.items() if k not in ("six", "weights")}
                        for r in records
                    ]
                ),
                mapModel=model,
                generator=generator,
                qualityGateVersion=VERSION,
                qualityGateSHA256=sha(gate_path),
            ),
            long_field_names=True,
            do_compression=True,
        )
    (folder / "preparation.json").write_text(
        json.dumps(
            dict(
                subject=subject,
                trainingCycles=9,
                independentTestCycles=0,
                allConditionsAreCalibration=True,
                MAPParametersRefitted=False,
                conditions=[
                    {
                        k: a[k]
                        for k in (
                            "id",
                            "trialPath",
                            "loadKg",
                            "speedCode",
                            "referenceValleys",
                            "torqueReferenceCoverageFraction",
                        )
                    }
                    for _, a in selected
                ],
                qualityGateVersion=VERSION,
                qualityGateSHA256=sha(gate_path),
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subjects", nargs="*", default=list(SUBJECTS))
    parser.add_argument("--output", type=Path, default=BASE)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--contact-at-peaks", action="store_true")
    args = parser.parse_args()
    for name in args.subjects:
        prepare(name, args.output, args.audit_only, args.contact_at_peaks)
