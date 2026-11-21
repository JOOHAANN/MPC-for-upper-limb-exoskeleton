"""Inspect all five calibration loads in new Exp folders without fitting models.

The historical joint-model checks are kept separate from file availability and
activation-only screening. Neither screen proves independent hardware sync.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import binary_closing, label

from audit_subjects import DATA, REQUIRED, ROOT, read_sto
from prepare_subjects import cycles
from quality_gate import (
    gaussian_matlab,
    gap_stats,
    native_tables,
    processed_signal,
    smoothness,
)

LOADS = (5, 7.5, 10, 12.5, 15)
SUFFIXES = (25, 35, 45, 55, 65)
PATTERN = re.compile(r"dynamics_stoop_(\d+(?:\.\d+)?)kg_(.+)")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def full_emg_screen(names, emg, ik):
    """Find candidate EMG bouts independently of the shorter CEINMS export.

    No activation values or missing samples are reconstructed. Interpolation is
    used to locate bouts; any bout near a >250 ms source gap is then rejected.
    Small signed sensor estimates are reported, not silently rectified for NMF.
    """
    issues, repairs = [], []
    if not np.isfinite(emg).all() or np.any(np.diff(emg[:, 0]) < 0):
        return dict(issues=["nonfinite_or_backward_filtered_emg"], bouts=[])
    times, first, counts = np.unique(emg[:, 0], return_index=True, return_counts=True)
    for i in np.flatnonzero(counts > 1):
        values = emg[emg[:, 0] == times[i], 1:]
        if not np.allclose(values, values[0], atol=1e-10, rtol=1e-9):
            issues.append("conflicting_duplicate_filtered_emg_timestamp")
        else:
            repairs.append("remove_identical_filtered_emg_duplicate_in_memory")
    if issues:
        return dict(issues=issues, repairs=repairs, bouts=[])
    source = emg[np.sort(first)]
    channels = [
        prefix + side
        for prefix in ("iliocost", "longlumb", "longthor")
        for side in ("_r", "_l")
    ]
    if not all(n in names for n in channels):
        return dict(issues=["missing_back_sensor_channels"], bouts=[])
    lo, hi = max(source[0, 0], ik[0, 0]), min(source[-1, 0], ik[-1, 0])
    if hi - lo < 1:
        return dict(issues=["insufficient_filtered_emg_ik_overlap"], bouts=[])
    fine = np.arange(0, hi - lo, 0.001)
    values = np.column_stack(
        [
            np.interp(fine, source[:, 0] - lo, source[:, names.index(n)])
            for n in channels
        ]
    )
    six = gaussian_matlab(gaussian_matlab(values, 5), 300)[::10]
    t, back = fine[::10], six.mean(axis=1)
    low, high = np.percentile(back, [10, 98])
    if high - low < 0.07:
        return dict(issues=["weak_full_record_back_emg_excursion"], bouts=[])
    active = back > low + 0.18 * (high - low)
    regions, count = label(binary_closing(active, structure=np.ones(71)))
    bouts = []
    for k in range(1, count + 1):
        ids = np.flatnonzero(regions == k)
        if len(ids) < 90:
            continue
        a, b = max(0, ids[0] - 25), min(len(t) - 1, ids[-1] + 25)
        if b - a > 1600:
            continue
        start, end = lo + t[a], lo + t[b]
        reasons = []
        maximum_gap = 0.0
        for clock in (source[:, 0], ik[:, 0]):
            delta = np.diff(clock)
            selected = delta[(clock[:-1] < end + 0.15) & (clock[1:] > start - 0.15)]
            maximum_gap = max(maximum_gap, float(max(selected, default=0)))
        if maximum_gap > 0.25:
            reasons.append("unsupported_emg_or_ik_gap_over_250ms")
        if active[0] and a < 70:
            reasons.append("left_censored_emg_bout")
        if active[-1] and b >= len(t) - 70:
            reasons.append("right_censored_emg_bout")
        part = six[a : b + 1]
        if part.min() < -0.05 or part.max() > 1.05:
            reasons.append("large_sensor_range_error")
        bouts.append(
            dict(
                startSecondsFromEmg=float(start - source[0, 0]),
                endSecondsFromEmg=float(end - source[0, 0]),
                durationSeconds=float(end - start),
                candidate=not reasons,
                concerns=reasons,
                maximumSourceGapSeconds=maximum_gap,
                minimumSmoothedSensor=float(part.min()),
                maximumSmoothedSensor=float(part.max()),
                smoothnessDiagnostic=smoothness(back[a : b + 1]),
            )
        )
    return dict(
        issues=[],
        repairs=repairs,
        commonEmgIkDurationSeconds=float(hi - lo),
        channelOrder=channels,
        bouts=bouts,
    )


def audit_trial(subject, path, limits):
    match = PATTERN.fullmatch(path.name)
    row = dict(
        subject=subject,
        trial=str(path.relative_to(DATA / subject)),
        loadKg=float(match[1]),
        speedCode=int(match[2]),
        status="audit_only",
    )
    result = dict(**row, nativeIssues=[], actions=[])
    names = (*REQUIRED, "emg.sto", "emgFilt.sto")
    result["missingFiles"] = [n for n in names if not (path / n).is_file()]
    result["fileSHA256"] = {n: sha(path / n) for n in names if (path / n).is_file()}
    result["sourceDigest"] = hashlib.sha256(
        json.dumps(result["fileSHA256"], sort_keys=True).encode()
    ).hexdigest()
    if result["missingFiles"]:
        return result, None
    try:
        tables, clocks, identity, exact, issues, repairs = native_tables(path)
        result.update(
            nativeIssues=issues,
            clocks=clocks,
            torqueIdentity=identity,
            activationEmgExactFraction=exact,
            repairs=repairs,
        )
        channels, emg = read_sto(path / "emgFilt.sto")
        delta = np.diff(emg[:, 0])
        result["filteredEmg"] = dict(
            channels=channels[1:],
            samples=len(emg),
            allFinite=bool(np.isfinite(emg).all()),
            backwardIntervals=int(sum(delta < 0)),
            duplicateTimestamps=int(len(emg) - len(np.unique(emg[:, 0]))),
            maximumGapSeconds=float(delta.max()),
            gapsOver250ms=int(sum(delta > 0.25)),
            durationSeconds=float(emg[-1, 0] - emg[0, 0]),
            startEpochSeconds=float(emg[0, 0]),
            endEpochSeconds=float(emg[-1, 0]),
            minimum=float(emg[:, 1:].min()),
            maximum=float(emg[:, 1:].max()),
        )
        ac = clocks["Activations.sto"]
        result["activationDurationSeconds"] = (
            ac["endEpochSeconds"] - ac["startEpochSeconds"]
        )
        result["activationEndBeforeFilteredEmgSeconds"] = float(
            emg[-1, 0] - ac["endEpochSeconds"]
        )
        result["fullEmgScreen"] = full_emg_screen(channels, emg, tables["ik.sto"][1])
        signal = processed_signal(path)
        result["commonDurationSeconds"] = float(signal["t"][-1])
        result["rawAngleRange"] = float(np.ptp(signal["qRaw"]))
        for cycle in cycles(row, signal):
            start = signal["absoluteStart"] + cycle["sourceStartSeconds"]
            end = signal["absoluteStart"] + cycle["sourceEndSeconds"]
            maximum_gap, gaps = gap_stats(tables, start - 0.15, end + 0.15)
            post = smoothness(cycle["back"])
            # These are activation-source checks inside the existing processing
            # chain, not a claim that every raw-EMG segment has been admitted.
            activation_reasons = list(issues)
            if cycle["six"].min() < -0.05 or cycle["six"].max() > 1.05:
                activation_reasons.append("large_channel_activation_range_error")
            if gaps:
                activation_reasons.append("unsupported_gap_over_250ms")
            if cycle["rightCensored"]:
                activation_reasons.append("incomplete_action_at_recording_end")
            rough = [key for key, limit in limits.items() if post[key] > limit]
            if rough:
                activation_reasons.append(
                    "activation_outside_benchmark_envelope:" + ",".join(rough)
                )
            joint_reasons = []
            if (
                min(cycle["back"].min(), cycle["front"].min()) < -1e-6
                or max(cycle["back"].max(), cycle["front"].max()) > 1 + 1e-6
            ):
                joint_reasons.append("lumped_hill_input_outside_zero_one")
            if cycle["angleTorqueCorrelation"] >= -0.25:
                joint_reasons.append("theta_torque_polarity_or_alignment")
            if np.ptp(cycle["q"]) < 0.2:
                joint_reasons.append("insufficient_theta_excursion")
            if result["rawAngleRange"] < 0.015:
                joint_reasons.append("raw_angle_range_too_small_for_minmax_scale")
            result["actions"].append(
                dict(
                    id=cycle["id"],
                    startSeconds=cycle["sourceStartSeconds"],
                    endSeconds=cycle["sourceEndSeconds"],
                    durationSeconds=cycle["duration"],
                    activationCandidate=not activation_reasons,
                    coupledCandidate=not (activation_reasons + joint_reasons),
                    activationConcerns=activation_reasons,
                    coupledConcerns=joint_reasons,
                    angleTorqueCorrelation=cycle["angleTorqueCorrelation"],
                    referenceValleyCount=len(cycle["referenceValleys"]),
                    maximumGapSeconds=maximum_gap,
                    gaps=gaps,
                    smoothnessAfter=post,
                    minimumChannelEstimate=float(cycle["six"].min()),
                    maximumChannelEstimate=float(cycle["six"].max()),
                )
            )
        return result, (channels, emg, signal)
    except Exception as exc:
        result["error"] = type(exc).__name__ + ": " + str(exc)
        return result, None


def save_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_examples(subject, examples, output):
    fig, axes = plt.subplots(5, 1, figsize=(12, 10), sharey=True)
    for ax, load in zip(axes, LOADS):
        row, (names, emg, signal) = examples[(subject, load)]
        t = emg[:, 0] - emg[0, 0]
        for prefix, label, color in (
            ("iliocost", "Iliocostalis", "#2166ac"),
            ("longlumb", "Lumbar longissimus", "#1b9e77"),
            ("longthor", "Thoracic longissimus", "#8c6bb1"),
        ):
            if all(prefix + side in names for side in ("_r", "_l")):
                y = (
                    emg[:, names.index(prefix + "_r")]
                    + emg[:, names.index(prefix + "_l")]
                ) / 2
                # Break the plotted line across acquisition gaps; do not draw an
                # interpolated bridge that could look like valid measurements.
                y[np.r_[False, np.diff(t) > 0.25]] = np.nan
                ax.plot(t[::8], y[::8], color=color, lw=0.8, label=label)
        ac = row["clocks"]["Activations.sto"]
        lo = ac["startEpochSeconds"] - emg[0, 0]
        hi = ac["endEpochSeconds"] - emg[0, 0]
        ax.axvspan(
            lo, hi, color=".5", alpha=0.12, label="Available activation interval"
        )
        ax.axvline(hi, color="#c34a36", ls="--", label="Activation export ends")
        ax.set(title=f"{load:g} kg | folder suffix 45", xlim=(0, t[-1]))
        ax.set_ylabel("Filtered EMG")
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("Native time from filtered EMG start (s)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(
        subject + " | all five loads: full EMG versus activation-export coverage"
    )
    fig.tight_layout(rect=(0, 0.065, 1, 0.97))
    fig.savefig(output / (subject + "_coverage.png"), dpi=145)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", nargs="+", default=["240725Exp", "240730Exp"])
    parser.add_argument(
        "--output", type=Path, default=ROOT / "docs/exp_data_audit_20261020"
    )
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    policy_file = ROOT / "docs/subject_quality_v2/quality_gate.json"
    limits = json.loads(policy_file.read_text())["smoothnessUpperLimits"]
    trials, inventory, summary, examples = [], [], [], {}
    for subject in args.subjects:
        folder = DATA / subject / "CalibrationData"
        for path in sorted(folder.glob("dynamics_stoop_*kg_*")):
            match = PATTERN.fullmatch(path.name)
            if not match or float(match[1]) not in LOADS:
                continue
            inventory.append(
                dict(
                    subject=subject,
                    trial=str(path.relative_to(DATA / subject)),
                    loadKg=float(match[1]),
                    suffix=match[2],
                    basicGrid=match[2] in {str(v) for v in SUFFIXES},
                )
            )
        for load in LOADS:
            selected = [
                r
                for r in inventory
                if r["subject"] == subject and r["loadKg"] == load and r["basicGrid"]
            ]
            for entry in selected:
                result, signals = audit_trial(
                    subject, DATA / subject / entry["trial"], limits
                )
                trials.append(result)
                if result["speedCode"] == 45 and signals is not None:
                    examples[(subject, load)] = result, signals
                print(
                    subject,
                    Path(entry["trial"]).name,
                    "native",
                    result.get("nativeIssues"),
                    "actions",
                    len(result["actions"]),
                    "activation",
                    sum(c["activationCandidate"] for c in result["actions"]),
                    "coupled",
                    sum(c["coupledCandidate"] for c in result["actions"]),
                    result.get("error", ""),
                    flush=True,
                )
            checked = [
                r for r in trials if r["subject"] == subject and r["loadKg"] == load
            ]
            actions = [c for r in checked for c in r["actions"]]
            full_bouts = [
                b for r in checked for b in r.get("fullEmgScreen", {}).get("bouts", [])
            ]
            summary.append(
                dict(
                    subject=subject,
                    loadKg=load,
                    basicFolders=len(checked),
                    completeFileSets=sum(not r["missingFiles"] for r in checked),
                    nativeChecksPassed=sum(
                        not r.get("nativeIssues")
                        and "error" not in r
                        and not r["missingFiles"]
                        for r in checked
                    ),
                    detectedActions=len(actions),
                    activationCandidates=sum(c["activationCandidate"] for c in actions),
                    coupledCandidates=sum(c["coupledCandidate"] for c in actions),
                    fullEmgCandidateBouts=sum(b["candidate"] for b in full_bouts),
                    fullEmgCandidateFiles=sum(
                        any(
                            b["candidate"]
                            for b in r.get("fullEmgScreen", {}).get("bouts", [])
                        )
                        for r in checked
                    ),
                    activationConcerns=dict(
                        Counter(v for c in actions for v in c["activationConcerns"])
                    ),
                    coupledConcerns=dict(
                        Counter(v for c in actions for v in c["coupledConcerns"])
                    ),
                )
            )
        if all((subject, load) in examples for load in LOADS):
            plot_examples(subject, examples, output)
    duplicates = []
    seen = {}
    for row in trials:
        digest = row["sourceDigest"]
        if digest in seen:
            duplicates.append([seen[digest], [row["subject"], row["trial"]]])
        seen[digest] = [row["subject"], row["trial"]]
    payload = dict(
        auditDate="2026-10-20",
        subjects=args.subjects,
        requiredLoadsKg=LOADS,
        folderSuffixes=SUFFIXES,
        suffixesAreNotIndependentlyMeasuredCadence=True,
        policySource=str(policy_file.relative_to(ROOT)),
        policySHA256=sha(policy_file),
        smoothnessUpperLimits=limits,
        hardwareSynchronizationProven=False,
        modelsFitted=False,
        modelParametersUnchanged=True,
        rawFilesModified=False,
        activationCandidatesAreNotFullRawEmgAdmission=True,
        fullEmgBoutsAreScreeningCandidatesNotTrainedActivations=True,
        duplicateBasicRecords=duplicates,
        inventory=inventory,
        summary=summary,
        trials=trials,
    )
    (output / "audit.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    save_csv(
        output / "load_summary.csv",
        [
            {k: (json.dumps(v) if isinstance(v, dict) else v) for k, v in r.items()}
            for r in summary
        ],
    )
    flat = [
        dict(
            subject=r["subject"],
            trial=r["trial"],
            loadKg=r["loadKg"],
            speedCode=r["speedCode"],
            missingFiles=";".join(r["missingFiles"]),
            nativeIssues=";".join(r["nativeIssues"]),
            error=r.get("error", ""),
            activationDurationSeconds=r.get("activationDurationSeconds"),
            filteredEmgDurationSeconds=r.get("filteredEmg", {}).get("durationSeconds"),
            activationEndsBeforeFilteredEmgSeconds=r.get(
                "activationEndBeforeFilteredEmgSeconds"
            ),
            filteredEmgMaximumGapSeconds=r.get("filteredEmg", {}).get(
                "maximumGapSeconds"
            ),
            detectedActions=len(r["actions"]),
            activationCandidates=sum(c["activationCandidate"] for c in r["actions"]),
            coupledCandidates=sum(c["coupledCandidate"] for c in r["actions"]),
        )
        for r in trials
    ]
    save_csv(output / "trial_summary.csv", flat)
    bouts = [
        dict(
            subject=r["subject"],
            trial=r["trial"],
            loadKg=r["loadKg"],
            speedCode=r["speedCode"],
            startSecondsFromEmg=b["startSecondsFromEmg"],
            endSecondsFromEmg=b["endSecondsFromEmg"],
            durationSeconds=b["durationSeconds"],
            candidate=b["candidate"],
            maximumSourceGapSeconds=b["maximumSourceGapSeconds"],
            concerns=";".join(b["concerns"]),
        )
        for r in trials
        for b in r.get("fullEmgScreen", {}).get("bouts", [])
    ]
    if bouts:
        save_csv(output / "full_emg_bouts.csv", bouts)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
