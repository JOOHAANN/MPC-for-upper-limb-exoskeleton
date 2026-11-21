"""Check clocks and activation against 240528 BEFORE any subject calibration.

No raw files are changed. Clock agreement is tested on native timestamps,
including the identity Torque = sum(MuscleForces * moment arms). This verifies
the exported CEINMS pipeline, not independently measured hardware latency.
The original 300 ms Gaussian activation filter is applied in physical time;
neither reference angles nor signal timing/amplitude are fitted to a model.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.io import loadmat
from scipy.signal import convolve, detrend, periodogram

from audit_subjects import DATA, OUT, ROOT, REQUIRED, read_sto, extract
from prepare_subjects import cycles

VERSION = "clock-smoothing-gate-v2"
MAX_GAP = 0.25


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nearest(times, reference):
    i = np.searchsorted(reference, times)
    return np.minimum(
        abs(times - reference[np.minimum(i, len(reference) - 1)]),
        abs(times - reference[np.maximum(i - 1, 0)]),
    )


def gaussian_matlab(y, window):
    """MATLAB's shrink-endpoint Gaussian, including even-window alignment."""
    # For an even window MATLAB uses [window/2, window/2-1] neighbours.
    # scipy's 'same' convolution with a symmetric even kernel does the same.
    x = np.arange(window) - (window - 1) // 2
    w = np.exp(-0.5 * (x / (window / 5)) ** 2)
    if y.ndim == 1:
        return convolve(y, w, mode="same") / convolve(np.ones(len(y)), w, mode="same")
    return np.column_stack([gaussian_matlab(v, window) for v in y.T])


def smoothness(y, dt=0.01):
    amplitude = max(float(np.percentile(y, 98) - np.percentile(y, 2)), 0.01)
    normalized = (y - np.mean(y)) / amplitude
    velocity = np.diff(normalized) / dt
    frequency, power = periodogram(detrend(normalized), fs=1 / dt, window="hann")
    return {
        "highFrequencyFraction": float(
            power[frequency >= 3].sum() / max(power.sum(), 1e-12)
        ),
        "normalizedSlopeP99": float(np.percentile(abs(velocity), 99)),
        "phaseTotalVariation": float(np.sum(abs(np.diff(normalized)))),
        "normalizedAccelerationRMS": float(
            np.sqrt(np.mean((np.diff(velocity) / dt) ** 2))
        ),
        "peak": float(max(y)),
    }


def native_tables(path):
    tables, clocks, issues, repairs = {}, {}, [], []
    for name in (*REQUIRED, "emg.sto"):
        names, data = read_sto(path / name)
        if len(data) < 30 or not np.isfinite(data).all():
            raise ValueError("Missing, short or non-finite data: " + name)
        delta = np.diff(data[:, 0])
        if np.any(delta < 0):
            issues.append("backward_clock:" + name)
        unique, first, count = np.unique(
            data[:, 0], return_index=True, return_counts=True
        )
        conflicts = 0
        for idx in np.flatnonzero(count > 1):
            values = data[data[:, 0] == unique[idx], 1:]
            conflicts += int(not np.allclose(values, values[0], atol=1e-10, rtol=1e-9))
        if conflicts:
            issues.append("conflicting_duplicate_timestamps:" + name)
        elif np.any(count > 1):
            repairs.append("remove_identical_duplicate_rows:" + name)
        clean = data[np.sort(first)]
        tables[name] = names, clean
        clocks[name] = {
            "samples": len(data),
            "nonpositiveIntervals": int(sum(delta <= 0)),
            "duplicateConflicts": conflicts,
            "medianDtSeconds": float(np.median(delta)),
            "maximumGapSeconds": float(max(delta)),
            "startEpochSeconds": float(data[0, 0]),
            "endEpochSeconds": float(data[-1, 0]),
        }
    anchor = tables["lmt.sto"][1][:, 0]
    for name, (_, table) in tables.items():
        supported = anchor[(anchor >= table[0, 0]) & (anchor <= table[-1, 0])]
        errors = nearest(supported, table[:, 0])
        clocks[name].update(
            exactAnchorFraction=float(np.mean(errors < 1e-6)),
            anchorNearestP95Seconds=float(np.percentile(errors, 95)),
            overlapFraction=float(len(supported) / len(anchor)),
        )
    # Native exact-time identity catches shifted/differently indexed exports.
    names, forces = tables["MuscleForces.sto"]
    arm_names, arms = tables["ma_L5_S1_Flex_Ext.sto"]
    _, torque = tables["Torque.sto"]
    common = np.intersect1d(np.intersect1d(forces[:, 0], arms[:, 0]), torque[:, 0])
    indices = [np.searchsorted(v[:, 0], common) for v in (forces, arms, torque)]
    matching = [arm_names.index(n) for n in names[1:]]
    rebuilt = np.sum(forces[indices[0], 1:] * arms[indices[1]][:, matching], axis=1)
    error = rebuilt - torque[indices[2], 1]
    identity = dict(
        samples=len(common),
        rmseNm=float(np.sqrt(np.mean(error**2))),
        maximumErrorNm=float(max(abs(error))),
        exactTimeFraction=len(common) / len(anchor),
    )
    if identity["rmseNm"] > 1e-4 or identity["exactTimeFraction"] < 0.98:
        issues.append("force_arm_torque_identity_failure")
    for name in ("ik.sto", "Torque.sto", "MuscleForces.sto", "ma_L5_S1_Flex_Ext.sto"):
        if clocks[name]["exactAnchorFraction"] < 0.98:
            issues.append("skeletal_clock_mismatch:" + name)
    activation, emg = tables["Activations.sto"][1], tables["emg.sto"][1]
    # Report small signed export/filter undershoots. Admission is based on
    # the actual lumped drive and the action's channels, not an unused idle
    # sample elsewhere in a long recording. Never silently rectify the data.
    clocks["Activations.sto"]["outOfRangeFraction"] = float(
        np.mean((activation[:, 1:] < -1e-6) | (activation[:, 1:] > 1.0 + 1e-6))
    )
    clocks["Activations.sto"]["minimumExportedActivation"] = float(
        activation[:, 1:].min()
    )
    clocks["Activations.sto"]["maximumExportedActivation"] = float(
        activation[:, 1:].max()
    )
    exact = float(np.mean(nearest(activation[:, 0], emg[:, 0]) < 1e-6))
    if exact < 0.99:
        issues.append("activation_emg_clock_mismatch")
    if clocks["Activations.sto"]["anchorNearestP95Seconds"] > 0.05:
        issues.append("activation_temporal_support_poor")
    return tables, clocks, identity, exact, issues, repairs


def processed_signal(path):
    original = extract(path)
    fine = np.arange(0, original["t"][-1] + 1e-7, 0.001)
    sample = np.arange(0, len(fine), 10)
    out = {"t": fine[sample], "absoluteStart": original["absoluteStart"]}
    for key in ("q", "torque", "qRaw"):
        # Keep the previously specified original Scale. No extra theta filter.
        out[key] = np.interp(out["t"], original["t"], original[key])
    for key in ("six", "back"):
        values = original[key]
        y = (
            np.interp(fine, original["t"], values)
            if values.ndim == 1
            else np.column_stack([np.interp(fine, original["t"], v) for v in values.T])
        )
        out[key + "Unfiltered"] = y[sample]
        # Same physical filter as the successful handwritten project.
        out[key] = gaussian_matlab(gaussian_matlab(y, 5), 300)[sample]
    out["front"] = np.interp(out["t"], original["t"], original["front"])
    out["weights"] = np.column_stack(
        [np.interp(out["t"], original["t"], v) for v in original["weights"].T]
    )
    return out


def gap_stats(tables, start, end):
    maximum, missing = 0.0, []
    for name in REQUIRED:
        t = tables[name][1][:, 0]
        gaps = np.diff(t)
        intersects = (t[:-1] < end) & (t[1:] > start)
        selected = gaps[intersects]
        maximum = max(maximum, float(max(selected, default=0)))
        bad = np.flatnonzero(intersects & (gaps > MAX_GAP))
        for i in bad:
            missing.append((max(start, float(t[i])), min(end, float(t[i + 1])), name))
    return maximum, missing


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            [
                {
                    k: (json.dumps(v) if isinstance(v, (dict, list)) else v)
                    for k, v in row.items()
                }
                for row in rows
            ]
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inventory",
        type=Path,
        default=ROOT / "results/multi_subject/data_quality.json",
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "results/multi_subject_quality_v2"
    )
    parser.add_argument("--loads", nargs="+", type=float, default=[5, 10, 15])
    parser.add_argument("--subjects", nargs="+")
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    (out / "signals").mkdir(exist_ok=True)
    baseline = []
    paired = loadmat(
        ROOT / "data/processed/hill_compliant/paired_trials_original_reference.mat",
        simplify_cells=True,
    )["trials"]
    for i, trial in enumerate(np.asarray(paired, dtype=object).reshape(-1)):
        baseline.append(
            dict(
                condition=i + 1,
                loadKg=float(trial["loadKg"]),
                **smoothness(trial["a"][:, 0]),
            )
        )
    # A broad benchmark envelope avoids treating genuine speed/amplitude
    # differences as noise. It is fixed before viewing new model outcomes.
    limits = {
        key: 2 * max(b[key] for b in baseline)
        for key in (
            "highFrequencyFraction",
            "normalizedSlopeP99",
            "phaseTotalVariation",
            "normalizedAccelerationRMS",
        )
    }
    inventory = json.loads(args.inventory.read_text(encoding="utf-8"))
    original_rows = [
        dict(
            subject="240528Pilot",
            trial=f"dynamics_stoop_{m}kg_{s}",
            loadKg=m,
            speedCode=s,
            status="baseline",
            sha256="baseline",
        )
        for m in (5, 10, 15)
        for s in (35, 45, 55)
    ]
    rows, decisions, plotted = [], [], []
    for row in original_rows + inventory:
        if row["loadKg"] not in args.loads:
            continue
        if args.subjects and row["subject"] not in args.subjects:
            continue
        result = {
            key: row[key]
            for key in ("subject", "trial", "loadKg", "speedCode", "status")
        }
        result.update(
            version=VERSION, issues=[], repairs=[], acceptedCycles=0, rejectedCycles=0
        )
        path = DATA / row["subject"] / row["trial"]
        result["fileSha256"] = {
            n: sha(path / n) for n in (*REQUIRED, "emg.sto") if (path / n).exists()
        }
        if row["status"] in ("excluded", "duplicate"):
            result["gate"] = "excluded"
            result["issues"] = [row["status"] + ":" + str(row.get("issues", []))]
            rows.append(result)
            continue
        try:
            tables, clocks, identity, ac_emg, issues, repairs = native_tables(path)
            result.update(
                clocks=clocks,
                torqueIdentity=identity,
                activationEmgExactFraction=ac_emg,
                issues=issues,
                repairs=repairs,
                appliedClockShiftSeconds=0.0,
            )
            s = processed_signal(path)
            candidates = cycles(row, s)
            for cycle in candidates:
                start = s["absoluteStart"] + cycle["sourceStartSeconds"]
                end = s["absoluteStart"] + cycle["sourceEndSeconds"]
                maxgap, gaps = gap_stats(tables, start - 0.15, end + 0.15)
                pre = smoothness(
                    s["backUnfiltered"][
                        (s["t"] >= cycle["sourceStartSeconds"])
                        & (s["t"] <= cycle["sourceEndSeconds"])
                    ]
                )
                post = smoothness(cycle["back"])
                reasons = list(issues)
                if (
                    min(cycle["back"].min(), cycle["front"].min()) < -1e-6
                    or max(cycle["back"].max(), cycle["front"].max()) > 1 + 1e-6
                ):
                    reasons.append("lumped_hill_input_outside_zero_one")
                if cycle["six"].min() < -0.05 or cycle["six"].max() > 1.05:
                    reasons.append("large_channel_activation_range_error")
                if gaps:
                    reasons.append("unsupported_gap_over_250ms")
                if cycle["rightCensored"]:
                    reasons.append("incomplete_action_at_recording_end")
                if cycle["angleTorqueCorrelation"] >= -0.25:
                    reasons.append("theta_torque_polarity_or_alignment")
                if np.ptp(cycle["q"]) < 0.2:
                    reasons.append("insufficient_theta_excursion")
                if np.ptp(s["qRaw"]) < 0.015:
                    reasons.append("raw_angle_range_too_small_for_minmax_scale")
                rough = [key for key, limit in limits.items() if post[key] > limit]
                if rough:
                    reasons.append(
                        "activation_outside_benchmark_envelope:" + ",".join(rough)
                    )
                decision = {
                    key: cycle[key]
                    for key in (
                        "id",
                        "trialPath",
                        "loadKg",
                        "speedCode",
                        "sourceStartSeconds",
                        "sourceEndSeconds",
                        "angleTorqueCorrelation",
                        "rightCensored",
                    )
                }
                decision.update(
                    subject=row["subject"],
                    accepted=not reasons,
                    reasons=reasons,
                    maximumGapSeconds=maxgap,
                    gaps=gaps,
                    smoothnessBefore=pre,
                    smoothnessAfter=post,
                    minimumFilteredChannelEstimate=float(cycle["six"].min()),
                    maximumFilteredChannelEstimate=float(cycle["six"].max()),
                    minimumActualHillSourceInput=float(
                        min(cycle["back"].min(), cycle["front"].min())
                    ),
                )
                decisions.append(decision)
                result["acceptedCycles"] += int(not reasons)
                result["rejectedCycles"] += int(bool(reasons))
            result["gate"] = "pass" if result["acceptedCycles"] else "excluded"
            result["signalCache"] = (
                row["subject"]
                + "_"
                + row["trial"].replace("\\", "_").replace("/", "_")
                + ".npz"
            )
            np.savez_compressed(out / "signals" / result["signalCache"], **s)
            result["signalCacheSHA256"] = sha(out / "signals" / result["signalCache"])
            plotted.append((result, s))
            if not candidates:
                result["issues"].append("no_complete_detectable_action")
        except Exception as exc:
            result["gate"] = "excluded"
            result["issues"].append(str(exc))
        rows.append(result)
        print(
            row["subject"],
            Path(row["trial"]).name,
            result["gate"],
            result["acceptedCycles"],
            result["rejectedCycles"],
            result["issues"],
            flush=True,
        )
    summary = []
    for subject in sorted(set(row["subject"] for row in rows)):
        subset = [d for d in decisions if d["subject"] == subject]
        passed = [d for d in subset if d["accepted"]]
        summary.append(
            dict(
                subject=subject,
                acceptedCycles=len(passed),
                rejectedCycles=len(subset) - len(passed),
                qualifiedLoadsKg=sorted(set(d["loadKg"] for d in passed)),
                rejectionReasons=dict(
                    Counter(reason for d in subset for reason in d["reasons"])
                ),
            )
        )
    payload = dict(
        version=VERSION,
        rawDataRoot=str(DATA),
        baselineSubject="240528Pilot",
        baseline=baseline,
        smoothnessUpperLimits=limits,
        maximumInterpolationGapSeconds=MAX_GAP,
        filterPolicy="Gaussian 5 ms then Gaussian 300 ms, sigma = window/5, centered OFFLINE",
        referencePolicy="Original full-record IK -> LMT -> Gaussian5 -> 0.5..1.4, unchanged",
        clockPolicy="Absolute timestamps, common support only, zero fitted clock shifts; exact-time pipeline identities",
        hardwareSynchronizationProven=False,
        evaluatedLoadsKg=args.loads,
        summary=summary,
        trials=rows,
        cycles=decisions,
    )
    (out / "quality_gate.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    write_csv(out / "cycle_decisions.csv", decisions)
    write_csv(out / "trial_clocks.csv", rows)
    with PdfPages(out / "data_quality_audit_v2.pdf") as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        subjects = [s["subject"] for s in summary]
        axes[0].bar(
            subjects,
            [s["acceptedCycles"] for s in summary],
            label="Qualified full action",
        )
        axes[0].bar(
            subjects,
            [s["rejectedCycles"] for s in summary],
            bottom=[s["acceptedCycles"] for s in summary],
            label="Excluded action",
            color="#dc8a61",
        )
        axes[0].tick_params(axis="x", rotation=25)
        axes[0].set_ylabel("Action segments")
        axes[0].legend(fontsize=8)
        groups = [d for d in decisions if d["subject"] != "240528Pilot"]
        axes[1].scatter(
            [d["smoothnessBefore"]["highFrequencyFraction"] * 100 for d in groups],
            [d["smoothnessAfter"]["highFrequencyFraction"] * 100 for d in groups],
            s=18,
            alpha=0.6,
        )
        axes[1].axhline(
            limits["highFrequencyFraction"] * 100,
            color="red",
            ls="--",
            label="Benchmark envelope",
        )
        axes[1].set(
            xlabel="Unfiltered activation: power above 3 Hz (%)",
            ylabel="Matched Gaussian activation: power above 3 Hz (%)",
        )
        axes[1].legend(fontsize=8)
        fig.suptitle(
            "Quality gate BEFORE subject calibration | no time shifts or angle rescaling"
        )
        fig.tight_layout()
        pdf.savefig(fig)
        fig.savefig(out / "quality_overview.png", dpi=160)
        plt.close(fig)
        for row, s in plotted:
            if row["subject"] == "240528Pilot":
                continue
            fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)
            axes[0].plot(
                s["t"],
                s["backUnfiltered"],
                color=".65",
                lw=0.7,
                label="Unfiltered weighted activation",
            )
            axes[0].plot(
                s["t"], s["back"], lw=1, label="Original 300 ms Gaussian policy"
            )
            axes[0].set_ylabel("Activation")
            axes[0].legend(fontsize=8)
            axes[1].plot(s["t"], s["q"], label="Original Scale reference")
            axes[1].set_ylabel("Theta (rad)")
            axes[2].plot(s["t"], s["torque"])
            axes[2].set(
                ylabel="Reference torque (N m)", xlabel="Common recording time (s)"
            )
            for d in decisions:
                if d["subject"] == row["subject"] and d["trialPath"] == row["trial"]:
                    for ax in axes:
                        ax.axvspan(
                            d["sourceStartSeconds"],
                            d["sourceEndSeconds"],
                            color="green" if d["accepted"] else "red",
                            alpha=0.08,
                        )
            for ax in axes:
                ax.grid(alpha=0.2)
            fig.suptitle(
                row["subject"]
                + " | "
                + row["trial"]
                + " | green = qualified, red = excluded",
                fontsize=10,
            )
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
