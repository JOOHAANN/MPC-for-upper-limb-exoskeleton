"""Select one observed M-shaped bout per condition and audit angle provenance.

This is a review, not model training. A later truncated repetition does not
invalidate an earlier complete bout. The user confirmed that these Exp records
have no Qualisys acquisition and authorized continued use of their existing IK.
"""

from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from audit_subjects import DATA, REQUIRED, ROOT, read_sto
from quality_gate import gaussian_matlab, native_tables, processed_signal, smoothness

SUBJECTS = ("240725Exp", "240730Exp")
LOADS = (5, 7.5, 10, 12.5, 15)
SPEEDS = (25, 35, 45, 55, 65)
OUTPUT = ROOT / "docs/exp_one_cycle_review"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def find_m_bouts(
    t, y, six, clocks, limits, check_signal=None, issues=(), excluded_intervals=()
):
    """Locate two main peaks in one bounded activity, without angle fitting."""
    baseline, high = np.percentile(y, [10, 98])
    span = high - baseline
    if span < 0.07:
        return []
    peaks, _ = find_peaks(
        y,
        height=baseline + 0.25 * span,
        prominence=max(0.025, 0.12 * span),
        distance=50,
        width=10,
    )
    candidates, seen = [], set()
    for a in range(len(peaks)):
        for b in range(a + 1, len(peaks)):
            p, q = int(peaks[a]), int(peaks[b])
            if not 0.6 <= t[q] - t[p] <= 10:
                continue
            amplitude = min(y[p], y[q]) - baseline
            level = baseline + 0.18 * amplitude
            valley = float(y[p : q + 1].min())
            # The inner valley must be distinct, without returning all the way
            # to the outer resting level between independent activities.
            if valley <= level or min(y[p], y[q]) - valley < max(
                0.025, 0.12 * amplitude
            ):
                continue
            left = np.flatnonzero(y[:p] <= level)
            right = np.flatnonzero(y[q + 1 :] <= level)
            if not len(left) or not len(right):
                continue
            lo, hi = int(left[-1]), int(q + 1 + right[0])
            if not 2 <= t[hi] - t[lo] <= 20:
                continue
            inside = peaks[(peaks >= lo) & (peaks <= hi)]
            # Ignore close shoulders on one lobe. Peak timing need not split
            # the window exactly in half: a real lift may have a long return.
            separated = [
                (int(u), int(v))
                for u in inside
                for v in inside
                if u < v
                and (t[v] - t[u]) >= 0.25 * (t[hi] - t[lo])
                and y[u : v + 1].min() > level
                and min(y[u], y[v]) - y[u : v + 1].min() >= max(0.025, 0.12 * amplitude)
            ]
            if not separated:
                continue
            strongest = list(max(separated, key=lambda pair: y[pair[0]] + y[pair[1]]))
            if (
                strongest != [p, q]
                or (lo, hi) in seen
                or (t[q] - t[p]) < 0.25 * (t[hi] - t[lo])
            ):
                continue
            seen.add((lo, hi))
            maximum_gap = 0.0
            for clock in clocks:
                gaps = np.diff(clock)
                selected = gaps[
                    (clock[:-1] < t[hi] + 0.15) & (clock[1:] > t[lo] - 0.15)
                ]
                maximum_gap = max(maximum_gap, float(max(selected, default=0)))
            concerns = list(issues)
            if any(
                t[lo] <= end and t[hi] >= start for start, end in excluded_intervals
            ):
                concerns.append("overlaps_conflicting_sensor_timestamp_exclusion")
            if maximum_gap > 0.25:
                concerns.append("source_gap_over_250ms_near_selected_bout")
            part = six[lo : hi + 1]
            if part.min() < -0.05 or part.max() > 1.05:
                concerns.append("large_channel_range_error")
            diagnostic = smoothness(
                (y if check_signal is None else check_signal)[lo : hi + 1]
            )
            rough = [key for key, limit in limits.items() if diagnostic[key] > limit]
            candidates.append(
                dict(
                    startEpochSeconds=float(t[lo]),
                    endEpochSeconds=float(t[hi]),
                    durationSeconds=float(t[hi] - t[lo]),
                    peakEpochSeconds=[float(t[p]), float(t[q])],
                    peakSeparationFraction=float((t[q] - t[p]) / (t[hi] - t[lo])),
                    peakValues=[float(y[p]), float(y[q])],
                    innerValley=valley,
                    otherSalientPeaks=int(len(inside) - 2),
                    maximumSourceGapSeconds=maximum_gap,
                    concerns=concerns,
                    smoothnessWarnings=rough,
                    smoothness=diagnostic,
                    cleanSourceWindow=not concerns,
                )
            )
    return candidates


def full_sensor_signal(path, ik, origin):
    names, data = read_sto(path / "emgFilt.sto")
    issues = []
    if not np.isfinite(data).all() or np.any(np.diff(data[:, 0]) < 0):
        return None, ["invalid_sensor_clock_or_values"]
    times, first, counts = np.unique(data[:, 0], return_index=True, return_counts=True)
    ambiguous = []
    for i in np.flatnonzero(counts > 1):
        rows = data[data[:, 0] == times[i], 1:]
        if not np.allclose(rows, rows[0], atol=1e-10, rtol=1e-9):
            ambiguous.append(float(times[i]))
            issues.append(
                "conflicting_duplicate_filtered_emg_timestamp_locally_excluded"
            )
    data = data[np.sort(first)]
    if ambiguous:
        # No conflicting values are averaged. Localization may interpolate the
        # removed point, but every accepted window excludes its filter support.
        data = data[~np.isin(data[:, 0], ambiguous)]
    channels = [
        prefix + side
        for prefix in ("iliocost", "longlumb", "longthor")
        for side in ("_r", "_l")
    ]
    lo, hi = max(data[0, 0], ik[0, 0]), min(data[-1, 0], ik[-1, 0])
    fine = np.arange(0, hi - lo, 0.001)
    values = np.column_stack(
        [np.interp(fine, data[:, 0] - lo, data[:, names.index(n)]) for n in channels]
    )
    six = gaussian_matlab(gaussian_matlab(values, 5), 300)[::10]
    return (
        dict(
            t=lo + fine[::10],
            six=six,
            y=six.mean(1),
            clock=data[:, 0],
            origin=origin,
            excludedIntervals=[(v - 0.2, v + 0.2) for v in ambiguous],
        ),
        issues,
    )


def provenance(subject):
    folder = DATA / subject
    configs = []
    for p in sorted(folder.glob("execution*.xml")):
        try:
            tree = ET.parse(p)
        except ET.ParseError:
            continue
        angles = [
            e.text.strip()
            for e in tree.iter()
            if e.tag.split("}")[-1] == "AngleDevice" and e.text
        ]
        if angles:
            configs.append(dict(path=str(p), sha256=sha(p), angleDevices=angles))
    matching = [
        str(p)
        for p in folder.rglob("*")
        if p.is_file()
        and (
            p.suffix.lower() in (".c3d", ".qtm", ".trc", ".mot")
            or (
                "qualisys" in p.name.lower()
                and p.suffix.lower() in (".mat", ".csv", ".txt")
            )
        )
    ]
    return dict(
        subject=subject,
        savedAngleDeviceConfigurations=configs,
        sameFolderOpticalSourceFiles=matching,
        currentIkStatus="existing_IK_accepted_by_user_Xsens_configured",
        runtimeSourceProven=False,
        qualisysReplacementRequired=False,
        userConfirmedNoQualisysForTheseExp=True,
        reason="The user confirmed no Qualisys data for these records and authorized their existing IK. Configuration points to Xsens; generic ik.sto headers alone do not identify the runtime device.",
    )


def review_trial(subject, path, limits):
    tables, clocks, identity, exact, issues, repairs = native_tables(path)
    _, raw = read_sto(path / "emgFilt.sto")
    origin = float(raw[0, 0])
    s = processed_signal(path)
    signal = dict(
        t=s["absoluteStart"] + s["t"], six=s["six"], y=s["six"].mean(1), origin=origin
    )
    activity = find_m_bouts(
        signal["t"],
        signal["y"],
        signal["six"],
        [a[:, 0] for _, a in tables.values()] + [raw[:, 0]],
        limits,
        check_signal=s["back"],
        issues=issues,
    )
    for r in activity:
        r["source"] = "CEINMS_activation"
    sensor, sensor_issues = full_sensor_signal(path, tables["ik.sto"][1], origin)
    raw_bouts = (
        []
        if sensor is None
        else find_m_bouts(
            sensor["t"],
            sensor["y"],
            sensor["six"],
            [sensor["clock"], tables["ik.sto"][1][:, 0]],
            limits,
            excluded_intervals=sensor["excludedIntervals"],
        )
    )
    for r in raw_bouts:
        r["source"] = "full_filtered_EMG_requires_activation_processing"
    all_bouts = activity + raw_bouts
    admissible = [r for r in all_bouts if r["cleanSourceWindow"]]
    # Selection is independent of fitted BSE/MPC outcomes. Prefer an existing
    # complete activation window, then a sensor-only restoration candidate.
    chosen = (
        min(
            admissible,
            key=lambda r: (
                r["source"] != "CEINMS_activation",
                bool(r["smoothnessWarnings"]),
                r["startEpochSeconds"],
            ),
        )
        if admissible
        else None
    )
    files = {n: sha(path / n) for n in (*REQUIRED, "emg.sto", "emgFilt.sto")}
    entry = dict(
        subject=subject,
        trial=str(path.relative_to(DATA / subject)),
        nativeIssues=issues,
        nativeRepairs=repairs,
        torqueIdentity=identity,
        activationEmgExactFraction=exact,
        sensorIssues=sensor_issues,
        sensorExcludedIntervals=[] if sensor is None else sensor["excludedIntervals"],
        fileSHA256=files,
        originEpochSeconds=origin,
        selected=chosen,
        candidates=all_bouts,
        modelsRetrained=False,
        kinematicsAcceptedByUser=True,
        qualisysVerified=False,
    )
    if chosen:
        chosen["startSecondsFromEmg"] = chosen["startEpochSeconds"] - origin
        chosen["endSecondsFromEmg"] = chosen["endEpochSeconds"] - origin
        chosen["peaksSecondsFromEmg"] = [t - origin for t in chosen["peakEpochSeconds"]]
        start, end = chosen["startEpochSeconds"], chosen["endEpochSeconds"]
        ik_names, ik = tables["ik.sto"]
        measured_angle = ik[
            (ik[:, 0] >= start) & (ik[:, 0] <= end), ik_names.index("L5_S1_Flex_Ext")
        ]
        chosen["nativeIkExcursion"] = (
            float(np.ptp(measured_angle)) if len(measured_angle) else 0.0
        )
        common_start = max(a[0, 0] for _, a in tables.values())
        common_end = min(a[-1, 0] for _, a in tables.values())
        chosen["existingMusculoskeletalTimeCoverageFraction"] = float(
            max(0.0, min(end, common_end) - max(start, common_start)) / (end - start)
        )
        lo, hi = max(start, common_start), min(end, common_end)
        chosen["angleTorqueCorrelationOnAvailableInterval"] = None
        if hi - lo > 1:
            index = (s["absoluteStart"] + s["t"] >= lo) & (
                s["absoluteStart"] + s["t"] <= hi
            )
            if (
                sum(index) > 10
                and np.std(s["q"][index]) > 1e-8
                and np.std(s["torque"][index]) > 1e-8
            ):
                chosen["angleTorqueCorrelationOnAvailableInterval"] = float(
                    np.corrcoef(s["q"][index], s["torque"][index])[0, 1]
                )
        index = (s["absoluteStart"] + s["t"] >= start) & (
            s["absoluteStart"] + s["t"] <= end
        )
        chosen["existingLumpedActivationInZeroOneOnAvailableInterval"] = (
            bool(
                s["back"][index].min() >= -1e-6
                and s["front"][index].min() >= -1e-6
                and s["back"][index].max() <= 1 + 1e-6
                and s["front"][index].max() <= 1 + 1e-6
            )
            if index.any()
            else None
        )
        display = signal if chosen["source"] == "CEINMS_activation" else sensor
    else:
        display = None
    return entry, display


def grid(subject, conditions, displays):
    fig, axes = plt.subplots(5, 5, figsize=(13, 11.5), sharey=True)
    for i, load in enumerate(LOADS):
        for j, speed in enumerate(SPEEDS):
            ax = axes[i, j]
            row = next(
                r
                for r in conditions
                if r["subject"] == subject
                and r["loadKg"] == load
                and r["speedCode"] == speed
            )
            choice = row["selected"]
            if choice:
                s = displays[(subject, load, speed)]
                ids = (s["t"] >= choice["startEpochSeconds"] - 0.15) & (
                    s["t"] <= choice["endEpochSeconds"] + 0.15
                )
                ax.plot(
                    s["t"][ids] - choice["startEpochSeconds"],
                    s["y"][ids],
                    color=(
                        "#1469b0"
                        if choice["source"] == "CEINMS_activation"
                        else "#c17a1b"
                    ),
                    lw=1.1,
                )
                ax.scatter(
                    np.asarray(choice["peakEpochSeconds"])
                    - choice["startEpochSeconds"],
                    choice["peakValues"],
                    s=13,
                    color="#993b32",
                    zorder=3,
                )
                tag = "A" if choice["source"] == "CEINMS_activation" else "E"
                suffix = ("*" if choice["smoothnessWarnings"] else "") + (
                    "R" if row.get("chosenFromNumberedRepeat") else ""
                )
                ax.set_title(f"{load:g} kg / {speed} | {tag}{suffix}", fontsize=9)
                ax.set_xlabel(
                    f"{choice['startSecondsFromEmg']:.1f}-{choice['endSecondsFromEmg']:.1f} s in source",
                    fontsize=7,
                )
            else:
                ax.set_facecolor("#eeeeee")
                ax.text(
                    0.5,
                    0.5,
                    "No clean M pair\nin this basic file",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=8,
                )
                ax.set_title(f"{load:g} kg / {speed}", fontsize=9)
            ax.grid(alpha=0.15)
    fig.suptitle(
        subject
        + " | one observed M-shaped pair per condition\nBlue A: activation; amber E: filtered EMG. Different processing stages; do not compare levels.",
        fontsize=12,
    )
    fig.text(
        0.5,
        0.015,
        "A: existing activation. E: full filtered EMG; activation still needs processing. *: smoothness review. R: numbered repeat. Dots: main peaks. No angle rescaling.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.945), h_pad=1.0)
    fig.savefig(OUTPUT / (subject + "_one_M_grid.png"), dpi=150)
    plt.close(fig)


def kept_example(conditions, displays):
    key = ("240725Exp", 5, 45)
    row = next(
        r for r in conditions if (r["subject"], r["loadKg"], r["speedCode"]) == key
    )
    chosen = row["selected"]
    if not chosen or chosen["source"] != "CEINMS_activation":
        return
    signal = displays[key]
    path = DATA / key[0] / row["trial"]
    names, emg = read_sto(path / "emgFilt.sto")
    channels = [
        prefix + side
        for prefix in ("iliocost", "longlumb", "longthor")
        for side in ("_r", "_l")
    ]
    measured = emg[:, [names.index(n) for n in channels]].mean(1)
    origin = row["originEpochSeconds"]
    lo, hi = chosen["startSecondsFromEmg"], chosen["endSecondsFromEmg"]
    end = signal["t"][-1] - origin
    fig, axes = plt.subplots(3, 1, figsize=(10, 7.2))
    axes[0].plot(emg[:, 0] - origin, measured, color="#c17a1b", lw=0.9)
    axes[1].plot(signal["t"] - origin, signal["y"], color="#3479ad", lw=1.1)
    for ax in axes[:2]:
        ax.axvspan(lo, hi, color="#39a16b", alpha=0.18)
        ax.axvline(end, color="#b33731", ls="--", lw=1.1)
        ax.set_xlim(0, 32)
    axes[1].axvspan(end, 32, color=".5", alpha=0.09)
    axes[1].text(
        24,
        0.5,
        "Later missing data\ndo not invalidate the kept M pair",
        transform=axes[1].get_xaxis_transform(),
        ha="center",
        fontsize=9,
    )
    keep = (signal["t"] >= chosen["startEpochSeconds"]) & (
        signal["t"] <= chosen["endEpochSeconds"]
    )
    axes[2].plot(signal["t"][keep] - origin, signal["y"][keep], color="#248553", lw=1.4)
    axes[2].scatter(
        chosen["peaksSecondsFromEmg"],
        chosen["peakValues"],
        color="#993b32",
        s=25,
        zorder=3,
    )
    axes[2].set_xlim(lo - 0.2, hi + 0.2)
    for ax, name in zip(
        axes,
        (
            "Full filtered EMG\nback-channel mean",
            "Available activation\nback-channel mean",
            "Retained activation\none M-shaped pair",
        ),
    ):
        ax.set_ylabel(name)
        ax.grid(alpha=0.18)
        ax.set_xlabel("Time from filtered-EMG start (s)")
    fig.suptitle(
        f"Keep the complete first M pair: {lo:.2f}-{hi:.2f} s\n240725Exp | 5 kg | folder suffix 45",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUTPUT / "240725Exp_keep_first_M.png", dpi=155)
    plt.close(fig)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    limits = json.loads(
        (ROOT / "docs/subject_quality_v2/quality_gate.json").read_text()
    )["smoothnessUpperLimits"]
    conditions, displays, source_audit, supplemental, basic_records = [], {}, [], [], []
    for subject in SUBJECTS:
        source_audit.append(provenance(subject))
        for load in LOADS:
            for speed in SPEEDS:
                path = (
                    DATA
                    / subject
                    / "CalibrationData"
                    / f"dynamics_stoop_{load:g}kg_{speed}"
                )
                row, display = review_trial(subject, path, limits)
                row.update(loadKg=load, speedCode=speed)
                conditions.append(row)
                basic_records.append(row)
                displays[(subject, load, speed)] = display
                c = row["selected"]
                print(
                    subject,
                    load,
                    speed,
                    c["source"] if c else "NO_CLEAN_PAIR",
                    (
                        (
                            round(c["startSecondsFromEmg"], 2),
                            round(c["endSecondsFromEmg"], 2),
                        )
                        if c
                        else ""
                    ),
                    flush=True,
                )
        for path in sorted(
            (DATA / subject / "CalibrationData").glob("dynamics_stoop_*kg_*_*")
        ):
            match = re.fullmatch(
                r"dynamics_stoop_(\d+(?:\.\d+)?)kg_(25|35|45|55|65)_(\d+)", path.name
            )
            if not match or float(match[1]) not in LOADS:
                continue
            try:
                row, repeat_display = review_trial(subject, path, limits)
            except Exception as exc:
                row = dict(
                    subject=subject,
                    trial=str(path.relative_to(DATA / subject)),
                    selected=None,
                    error=type(exc).__name__ + ": " + str(exc),
                )
            row.update(
                loadKg=float(match[1]),
                speedCode=int(match[2]),
                includedInBasicGrid=False,
                numberedRepeatLabelConvention=True,
            )
            supplemental.append(row)
            target = next(
                i
                for i, r in enumerate(conditions)
                if r["subject"] == subject
                and r["loadKg"] == row["loadKg"]
                and r["speedCode"] == row["speedCode"]
            )
            if not conditions[target]["selected"] and row["selected"]:
                conditions[target] = dict(
                    **row,
                    chosenFromNumberedRepeat=True,
                    replacedBasicTrial=conditions[target]["trial"],
                )
                displays[(subject, row["loadKg"], row["speedCode"])] = repeat_display
            print("SUPPLEMENTAL", subject, path.name, bool(row["selected"]), flush=True)
        grid(subject, conditions, displays)
    summary = []
    for subject in SUBJECTS:
        rows = [r for r in conditions if r["subject"] == subject]
        summary.append(
            dict(
                subject=subject,
                conditions=25,
                basicSelectedConditions=sum(
                    r["selected"] is not None
                    for r in basic_records
                    if r["subject"] == subject
                ),
                numberedRepeatRescues=sum(
                    bool(r.get("chosenFromNumberedRepeat")) for r in rows
                ),
                selectedConditions=sum(r["selected"] is not None for r in rows),
                sourceCounts=dict(
                    Counter(
                        r["selected"]["source"] if r["selected"] else "no_clean_pair"
                        for r in rows
                    )
                ),
                smoothnessReviewConditions=sum(
                    bool(r["selected"] and r["selected"]["smoothnessWarnings"])
                    for r in rows
                ),
                existingMusculoskeletalWindowComplete=sum(
                    bool(
                        r["selected"]
                        and r["selected"]["existingMusculoskeletalTimeCoverageFraction"]
                        > 0.999
                    )
                    for r in rows
                ),
                kinematicsAcceptedByUser=True,
                modelsRetrained=False,
                missing=[
                    dict(loadKg=r["loadKg"], speedCode=r["speedCode"])
                    for r in rows
                    if not r["selected"]
                ],
            )
        )
    payload = dict(
        rule="One complete observed M pair is sufficient for condition-level data inclusion; later partial bouts do not reject earlier complete bouts.",
        hardwareSynchronizationProven=False,
        qualisysRequired=False,
        userConfirmedNoQualisysForTheseExp=True,
        kinematicsPolicy="Use existing recorded IK; retain the original reference Scale when preparing a new model. No angle sign flip or selected-bout rescaling is performed here.",
        modelsRetrained=False,
        originalDataModified=False,
        rawSensorIsNotActivation=True,
        summary=summary,
        sourceAudit=source_audit,
        conditions=conditions,
        basicRecords=basic_records,
        supplementalNumberedRecords=supplemental,
    )
    (OUTPUT / "review.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    flat = []
    for r in conditions:
        c = r["selected"] or {}
        flat.append(
            dict(
                subject=r["subject"],
                loadKg=r["loadKg"],
                speedCode=r["speedCode"],
                trial=r["trial"],
                source=c.get("source", "no_clean_pair"),
                startSecondsFromEmg=c.get("startSecondsFromEmg"),
                endSecondsFromEmg=c.get("endSecondsFromEmg"),
                durationSeconds=c.get("durationSeconds"),
                maximumSourceGapSeconds=c.get("maximumSourceGapSeconds"),
                existingMusculoskeletalTimeCoverageFraction=c.get(
                    "existingMusculoskeletalTimeCoverageFraction"
                ),
                nativeIkExcursion=c.get("nativeIkExcursion"),
                angleTorqueCorrelationOnAvailableInterval=c.get(
                    "angleTorqueCorrelationOnAvailableInterval"
                ),
                smoothnessWarnings=";".join(c.get("smoothnessWarnings", [])),
                kinematicsAcceptedByUser=True,
                modelsRetrained=False,
            )
        )
    with (OUTPUT / "conditions.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    (OUTPUT / "kinematics_policy.json").write_text(
        json.dumps(
            dict(
                qualisysReplacementRequired=False,
                userConfirmedNoQualisysForTheseExp=True,
                currentIKAccepted=True,
                note="Continue with existing ik.sto per the user's correction. Keep device provenance explicit; do not substitute another Pilot session's optical data.",
                records=[
                    dict(
                        subject=r["subject"],
                        trial=r["trial"],
                        ikPath=str(DATA / r["subject"] / r["trial"] / "ik.sto"),
                        ikSHA256=r["fileSHA256"]["ik.sto"],
                    )
                    for r in conditions
                ],
            ),
            indent=2,
        )
        + "\n"
    )
    kept_example(conditions, displays)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
