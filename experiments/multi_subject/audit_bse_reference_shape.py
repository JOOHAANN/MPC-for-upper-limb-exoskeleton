"""Trace the Exp reference waveforms without altering data or fitted models.

This is an acquisition/legacy-code audit, not a torque recalibration. Outputs
are kept beside the MATLAB diagnostics. Gaussian-filtered torque is explicitly
display-only and never replaces a reference used for fitting or scoring.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import find_peaks

from audit_subjects import DATA, ROOT, read_sto, extract, clock_data, REPRESENTATIVES
from compare_exp_signal_levels import SENSOR_CHANNELS
from control_nine_cases import many
from prepare_exp_three_loads import angle_reference, supported
from recompute_exp_activation import selected_signal, MODEL_CHANNELS
from quality_gate import gaussian_matlab, sha

OUT = ROOT / "Demo/BSE/diagnostics"
LEGACY = ROOT.parent / "YOUHAN LI code and project description"
SUBJECTS = ("240725Exp_A01", "240730Exp_A01")


def correlation(a, b):
    keep = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[keep], b[keep])[0, 1]) if keep.sum() > 3 else None


def normalize(y):
    return y / max(float(np.nanmax(abs(y))), 1e-9)


def smooth_observed(y):
    result = np.full_like(y, np.nan)
    changes = np.diff(np.r_[False, np.isfinite(y), False].astype(int))
    for start, end in zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)):
        result[start:end] = gaussian_matlab(y[start:end], min(30, end - start))
    return result


def missing_spans(t, valid):
    changes = np.diff(np.r_[False, ~valid, False].astype(int))
    for start, end in zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)):
        yield t[start], t[min(end, len(t) - 1)]


def legacy_archives():
    fig, axes = plt.subplots(3, 2, figsize=(15, 9), constrained_layout=True)
    audits = []
    for column, tag in enumerate(("B", "C")):
        folder = LEGACY / f"My_MAP_model - 5 10 15 - subject {tag}"
        combination = loadmat(folder / "Combination.mat", simplify_cells=True)[
            "Combination"
        ]
        raw_path = folder / "Combination_raw.mat"
        raw = (
            loadmat(raw_path, simplify_cells=True)["Combination"]
            if raw_path.exists()
            else combination
        )
        t = np.arange(len(combination)) * 0.001
        axes[0, column].plot(
            t, combination[:, 0], color="#d04838", label="Saved torque reference"
        )
        for row, filename, variable in (
            (0, "Torque_Tac0.mat", "Torque_Tac0"),
            (2, "Theta_Tac0.mat", "Theta_Tac0"),
        ):
            values = np.asarray(
                loadmat(folder / filename, simplify_cells=True)[variable]
            ).reshape(-1)
            if len(values) == len(t):
                axes[row, column].plot(
                    t,
                    values,
                    color="#176bb0",
                    alpha=0.85,
                    label="Saved legacy model output",
                )
        axes[1, column].plot(
            t,
            raw[:, 1],
            color="#178561",
            label=(
                "Saved pre-MAP profile"
                if raw_path.exists()
                else "Saved prepared profile"
            ),
        )
        if raw_path.exists():
            axes[1, column].plot(
                t, combination[:, 1], color="#7950a4", label="Saved MAP profile"
            )
        axes[2, column].plot(
            t, combination[:, 2], color="#d04838", label="Saved angle reference"
        )
        for row, label in enumerate(
            ("Torque (N m)", "Prepared activation", "Theta (rad)")
        ):
            ax = axes[row, column]
            ax.set_ylabel(label)
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8, loc="upper right")
            ax.set_xlim(t[[0, -1]])
        axes[2, column].set_xlabel("Saved concatenated time (s)")
        axes[0, column].set_title(
            (
                "Legacy B | code reads 240725Exp"
                if tag == "B"
                else "Legacy C | code paths read 240528Pilot, not 240730Exp"
            ),
            fontsize=11,
        )
        audits.append(
            dict(
                folder=str(folder),
                codeSubject="240725Exp" if tag == "B" else "240528Pilot",
                savedRows=len(combination),
                savedReferenceContainsNaN=bool(np.isnan(combination[:, 0]).any()),
                archiveSHA256=sha(folder / "Combination.mat"),
                waveformSubjectProvenanceConfirmed=False,
            )
        )
        if tag == "C":
            geometry = loadmat(
                ROOT / "data/processed/demo_calibration__geometry/Combination.mat",
                simplify_cells=True,
            )["Combination"]
            audits[-1]["savedTorqueAndThetaMatchArchivedGeometryCombination"] = bool(
                combination.shape == geometry.shape
                and np.array_equal(combination[:, [0, 2]], geometry[:, [0, 2]])
            )
    fig.suptitle(
        "Archived curves are processed combinations, not untouched acquisition references",
        fontsize=15,
    )
    fig.savefig(OUT / "legacy_saved_curves.png", dpi=160)
    plt.close(fig)
    return audits


def old_b_windows():
    """Translate the handwritten concatenation indices back to sensor clocks."""
    offset = 0.0
    rows = []
    plots = []
    for load, endpoint, start, end in (
        (5, 500, 4361, 10796),
        (10, 492, 19954, 26663),
        (15, 435, 31924, 38728),
    ):
        path = DATA / "240725Exp" / "CalibrationData" / f"dynamics_stoop_{load}kg_45"
        _, lmt = read_sto(path / "lmt.sto")
        _, sensor = read_sto(path / "emgFilt.sto")
        first = lmt[0, 0]
        finish = lmt[endpoint - 1, 0]
        origin = sensor[0, 0]
        # The old code has one-millisecond sampling and strictly open windows.
        lo = first - origin + start * 0.001 - offset
        hi = first - origin + (end - 2) * 0.001 - offset
        count = int(np.floor((finish - first) / 0.001 + 0.0003)) + 1
        clock = origin + np.arange(lo, hi + 0.0001, 0.01)
        _, torque = read_sto(path / "Torque.sto")
        valid = supported(torque[:, 0], clock)
        angle_issue = ""
        try:
            q, _ = angle_reference(path, clock)
        except ValueError as error:
            # Show only historically supported samples when the handwritten
            # window fails today's continuity gate. Do not repair that gap.
            angle_issue = str(error)
            historical = extract(path)
            native = historical["absoluteStart"] + historical["t"]
            q = np.interp(clock, native, historical["q"])
            _, ik = read_sto(path / "ik.sto")
            q[~(supported(ik[:, 0], clock) & supported(native, clock))] = np.nan
        y = np.interp(clock, torque[:, 0], torque[:, 1])
        y[~valid] = np.nan
        row = dict(
            loadKg=load,
            speedCode=45,
            startSecondsFromSensor=lo,
            endSecondsFromSensor=hi,
            originalTorqueCoverageFraction=float(valid.mean()),
            angleContinuousSupportPassed=not bool(angle_issue),
            angleIssue=angle_issue,
            angleValleys=len(find_peaks(-q, prominence=0.12, distance=50)[0]),
            missingTorqueIntervalsSecondsFromSensor=[
                list(span) for span in missing_spans(clock - origin, valid)
            ],
        )
        try:
            item, _ = selected_signal(
                dict(
                    subject="240725Exp",
                    trial=str(path.relative_to(DATA / "240725Exp")),
                    originEpochSeconds=origin,
                    selected=dict(
                        source="filtered",
                        startEpochSeconds=clock[0],
                        endEpochSeconds=clock[-1],
                    ),
                ),
                {name: 0.1 for name in MODEL_CHANNELS},
            )
            row["excitationWindowPassedExistingSignalChecks"] = True
        except ValueError as error:
            row["excitationWindowPassedExistingSignalChecks"] = False
            row["inputIssue"] = str(error)
        rows.append(row)
        plots.append(dict(t=clock - clock[0], torque=y, q=q))
        offset += count * 0.001
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), constrained_layout=True)
    for k, (row, data) in enumerate(zip(rows, plots)):
        axes[0, k].plot(
            data["t"], data["torque"], color="#d04838", label="Original CEINMS torque"
        )
        axes[1, k].plot(
            data["t"], data["q"], color="#464f5d", label="Original IK scale"
        )
        for axis, signal in zip(axes[:, k], (data["torque"], data["q"])):
            for lo, hi in missing_spans(data["t"], np.isfinite(signal)):
                axis.axvspan(lo, hi, color=".7", alpha=0.25)
        axes[0, k].set_title(
            f"{row['loadKg']} kg / 45 | sensor {row['startSecondsFromSensor']:.2f}-{row['endSecondsFromSensor']:.2f} s\nTorque support {100*row['originalTorqueCoverageFraction']:.1f}%; IK continuity {'OK' if row['angleContinuousSupportPassed'] else 'fails'}"
        )
        for axis in axes[:, k]:
            axis.grid(alpha=0.2)
            axis.set_xlabel("Time inside legacy window (s)")
        axes[0, k].set_ylabel("Torque (N m)")
        axes[1, k].set_ylabel("Theta (rad)")
    fig.suptitle(
        "Handwritten B windows on the raw acquisition clock | no amplitude scaling or reference filling",
        fontsize=13,
    )
    fig.savefig(OUT / "legacy_b_raw_windows.png", dpi=160)
    plt.close(fig)
    return rows


def check_native_input_clocks():
    """Compare mapped raw-sensor channels with the native logged excitation."""
    rows = []
    for subject, trial in (
        ("240725Exp", "dynamics_stoop_5kg_45"),
        ("240730Exp", "dynamics_stoop_5kg_45_3"),
    ):
        path = DATA / subject / "CalibrationData" / trial
        sn, sensor = read_sto(path / "emgFilt.sto")
        en, logged = read_sto(path / "emg.sto")
        sensor = clock_data(sensor)
        logged = clock_data(logged)
        start = max(sensor[0, 0], logged[0, 0]) + 0.3
        finish = min(sensor[-1, 0], logged[-1, 0]) - 0.3
        t = start + np.arange(int((finish - start) / 0.001)) * 0.001
        valid = supported(logged[:, 0], t) & supported(sensor[:, 0], t)
        for sensor_name, muscle_name in zip(SENSOR_CHANNELS, REPRESENTATIVES):
            observed = gaussian_matlab(
                np.interp(t, logged[:, 0], logged[:, en.index(muscle_name)]), 300
            )

            def score(shift):
                raw = gaussian_matlab(
                    np.interp(
                        t + shift, sensor[:, 0], sensor[:, sn.index(sensor_name)]
                    ),
                    300,
                )
                return correlation(
                    np.where(valid, raw, np.nan), np.where(valid, observed, np.nan)
                )

            shifts = np.arange(-0.25, 0.2501, 0.01)
            scores = np.array([score(shift) for shift in shifts])
            rows.append(
                dict(
                    subject=subject,
                    trial=trial,
                    sensorChannel=sensor_name,
                    modelChannel=muscle_name,
                    sameTimeCorrelation=score(0),
                    bestSensorQueryOffsetSeconds=float(shifts[np.argmax(scores)]),
                    bestCorrelation=float(scores.max()),
                )
            )
    return rows


def audit_subject(subject, output):
    base = ROOT / "results/nine_case_75pct" / subject
    prepared = loadmat(base / "prepared.mat", simplify_cells=True)
    gate = json.loads((base / "quality_gate.json").read_text())
    cases = many(
        next(s for s in many(output["subjects"]) if s["id"] == subject)["cases"]
    )
    rows = []
    details = []
    fig, axes = plt.subplots(3, 3, figsize=(15, 10), constrained_layout=True)
    for index, (cycle, audit, case) in enumerate(
        zip(many(prepared["cycles"]), gate["trials"], cases)
    ):
        original = DATA / subject.removesuffix("_A01") / cycle["trialPath"]
        derived = DATA / subject / cycle["trialPath"]
        _, old = read_sto(original / "Torque.sto")
        _, new = read_sto(derived / "Torque.sto")
        _, sensor = read_sto(original / "emgFilt.sto")
        origin = sensor[0, 0]
        t = np.arange(len(cycle["q"])) * 0.01
        epoch = origin + cycle["sourceStartSeconds"] + t
        observed = np.asarray(cycle["torqueReferenceValid"], dtype=bool)
        torque = np.asarray(cycle["torque"])
        smoothed = smooth_observed(torque)
        older = np.interp(epoch, old[:, 0], old[:, 1])
        older[~supported(old[:, 0], epoch)] = np.nan
        a = np.asarray(cycle["back"])
        model = case["outputs"][:, 2, 0]
        common, oi, ni = np.intersect1d(old[:, 0], new[:, 0], return_indices=True)
        outside = (epoch < new[0, 0]) | (epoch > new[-1, 0])
        stage_correlations = {
            name: correlation(case["outputs"][:, column, 0], torque)
            for name, column in (
                ("backNeuralCommand", 17),
                ("effectiveBackActivation", 13),
                ("backTorque", 3),
                ("abdominalTorque", 4),
                ("netTorque", 2),
            )
        }
        row = dict(
            subject=subject,
            loadKg=float(cycle["loadKg"]),
            speedCode=int(cycle["speedCode"]),
            recordId=cycle["id"],
            inputSource=audit["source"],
            selectedStartSeconds=float(cycle["sourceStartSeconds"]),
            selectedEndSeconds=float(cycle["sourceEndSeconds"]),
            sensorLastSeconds=float(sensor[-1, 0] - origin),
            torqueFirstSeconds=float(new[0, 0] - origin),
            torqueLastSeconds=float(new[-1, 0] - origin),
            torqueNativeMedianDt=float(np.median(np.diff(new[:, 0]))),
            maximumTorqueGapSeconds=float(np.diff(new[:, 0]).max()),
            torqueCoverageFraction=float(observed.mean()),
            activationReferenceCorrelation=correlation(a, torque),
            outsideNativeTorqueFraction=float(outside.mean()),
            interiorTorqueGapFraction=float((~observed & ~outside).mean()),
            bseExperimentalReferenceCorrelation=correlation(model, torque),
            oldNewTorqueCorrelation=correlation(old[oi, 1], new[ni, 1]),
            referenceSHA256=sha(derived / "Torque.sto"),
        )
        rows.append(row)
        details.append(
            dict(
                time=t,
                activation=a,
                reference=torque,
                referenceDisplaySmoothed=smoothed,
                oldReference=older,
                backTorque=case["outputs"][:, 3, 0],
                abdomenTorque=case["outputs"][:, 4, 0],
                netTorque=model,
                loadKg=row["loadKg"],
                speedCode=row["speedCode"],
                neuralCommand=case["outputs"][:, 17, 0],
                effectiveHillActivation=case["outputs"][:, 13, 0],
                stageCorrelations=stage_correlations,
                missingTorqueIntervals=[
                    list(span) for span in missing_spans(t, observed)
                ],
            )
        )
        ax = axes.flat[index]
        ax.plot(
            t,
            normalize(a),
            color="#148865",
            lw=1.8,
            label="Actual experimental input, A=0.1",
        )
        ax.plot(
            t,
            normalize(torque),
            color="#d04838",
            lw=0.7,
            alpha=0.35,
            label="CEINMS estimate, unsmoothed",
        )
        ax.plot(
            t,
            normalize(smoothed),
            color="#d04838",
            lw=1.6,
            label="CEINMS, 300 ms display filter",
        )
        ax.plot(
            t,
            normalize(model),
            color="#176bb0",
            lw=1.4,
            label="BSE net torque: experimental input",
        )
        for lo, hi in missing_spans(t, observed):
            ax.axvspan(lo, hi, color=".7", alpha=0.25)
        ax.set_title(
            f"{row['loadKg']:g} kg / {row['speedCode']} | torque support {100*observed.mean():.0f}%\nr(input, ref)={row['activationReferenceCorrelation']:.3f}; r(BSE, ref)={row['bseExperimentalReferenceCorrelation']:.3f}",
            fontsize=10,
        )
        ax.set_xlabel("Time inside selected action (s)")
        ax.set_ylabel("Own-peak normalized amplitude")
        ax.grid(alpha=0.2)
        ax.set_xlim(t[[0, -1]])
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, fontsize=9)
    fig.suptitle(
        subject
        + " | the blue BSE torque is not the experimental activation signal\nGray = unavailable reference; normalization and smoothing are display-only.",
        fontsize=13,
    )
    short = subject[3:6]
    fig.savefig(OUT / f"{short}_input_reference_comparison.png", dpi=160)
    plt.close(fig)
    savemat(
        OUT / f"{short}_reference_audit.mat",
        dict(cases=np.array(details, dtype=object)),
        do_compression=True,
        long_field_names=True,
    )
    (OUT / f"{short}_signal_path_audit.json").write_text(
        json.dumps(
            [
                {
                    k: v
                    for k, v in item.items()
                    if k
                    in (
                        "loadKg",
                        "speedCode",
                        "stageCorrelations",
                        "missingTorqueIntervals",
                    )
                }
                for item in details
            ],
            indent=2,
        )
        + "\n"
    )
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = loadmat(ROOT / "Demo/BSE/output/bse_results.mat", simplify_cells=True)[
        "results"
    ]
    legacy = legacy_archives()
    windows = old_b_windows()
    clocks = check_native_input_clocks()
    rows = [row for subject in SUBJECTS for row in audit_subject(subject, results)]
    with (OUT / "reference_quality.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = dict(
        referenceCurveIsModelEstimate=True,
        existingModelsChanged=False,
        activationCoverageIsNotTorqueCoverage=True,
        referenceDashPatternWasPlotStyle=True,
        displaySmoothingChangesTrainingTargets=False,
        legacyArchives=legacy,
        legacyBRawWindows=windows,
        nativeExcitationClockChecks=clocks,
        conditions=rows,
    )
    (OUT / "reference_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                conditions=len(rows),
                fullTorqueConditions=sum(
                    r["torqueCoverageFraction"] == 1 for r in rows
                ),
                minimumInputReferenceCorrelation=min(
                    r["activationReferenceCorrelation"] for r in rows
                ),
                legacyBRawWindows=windows,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
