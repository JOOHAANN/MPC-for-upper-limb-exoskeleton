"""Train the fourth subject MAP and select nine measured BSE conditions.

MAP uses excitation only. BSE additionally requires measured IK and at least
some supported CEINMS torque. Raw files and the original angle scale stay intact.
Selection precedes BSE fitting and does not use a fitted model's performance.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import find_peaks

from audit_subjects import DATA, ROOT, REPRESENTATIVES, read_sto
from control_nine_cases import digest, many
from excitation_only_map import build_review
from fit_hill_compliant import input_generator
from prepare_exp_three_loads import angle_reference, supported
from prepare_subjects import combine, fit_map, map_predict, structs
from quality_gate import native_tables

SUBJECT = "240617Pilot"
OUT = ROOT / "results/four_subject_demo" / SUBJECT


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gate, review = build_review(SUBJECT, OUT / "excitation")
    old = loadmat(
        ROOT / "results/multi_subject_quality_v2" / SUBJECT / "prepared.mat",
        simplify_cells=True,
    )
    weights = np.asarray(old["mapModel"]["lumpWeights"])
    shape_audits = []
    records = []
    for row in review["conditions"]:
        path = DATA / SUBJECT / row["trial"]
        # Confirm that +0.1 is the existing processing, not a new label scale.
        names, emg = read_sto(path / "emg.sto")
        anames, activation = read_sto(path / "Activations.sto")
        common, ei, ai = np.intersect1d(
            emg[:, 0], activation[:, 0], return_indices=True
        )
        errors = []
        for channel in (*REPRESENTATIVES, "rect_abd_r", "rect_abd_l"):
            x = emg[ei, names.index(channel)]
            y = activation[ai, anames.index(channel)]
            errors.extend(np.expm1(0.1 * x) / np.expm1(0.1) - y)
        error = float(np.sqrt(np.mean(np.square(errors))))
        if len(common) < 30 or error > 1e-6:
            raise ValueError("Native activation shape differs from +0.1: " + str(path))
        shape_audits.append(
            dict(trial=row["trial"], samples=len(common), responseRMSE=error)
        )
        item = loadmat(row["windowPath"], simplify_cells=True)
        t = np.asarray(item["timeEpochSeconds"])
        six = np.asarray(item["activation"])[:, :6]
        duration = float(t[-1] - t[0])
        phase_time = np.linspace(0, duration, 200)
        profile = np.column_stack([np.interp(phase_time, t - t[0], v) for v in six.T])
        records.append(
            dict(
                id=f"{row['loadKg']:g}kg_{row['speedCode']}",
                trialPath=row["trial"],
                loadKg=row["loadKg"],
                speedCode=row["speedCode"],
                cadenceBPM=60 / duration,
                duration=duration,
                profile=profile,
                six=six,
                back=six @ weights,
                front=np.asarray(item["activation"])[:, 6],
                q=np.zeros(len(t)),
                weights=np.tile(weights, (len(t), 1)),
                epoch=t,
                split="train",
                sourceWindowPath=row["windowPath"],
            )
        )
    # Hold out an entire load/speed condition per load before MAP fitting.
    for load in (5, 7.5, 10, 12.5, 15):
        group = sorted(
            [r for r in records if r["loadKg"] == load], key=lambda r: r["speedCode"]
        )
        if len(group) < 2:
            raise ValueError(
                f"Insufficient excitation-only MAP coverage at {load:g} kg"
            )
        group[-1]["split"] = "test_condition"
    train = [r for r in records if r["split"] == "train"]
    model = fit_map(train)
    map_rows = []
    for r in records:
        r["mappedBack"] = map_predict(model, r)
        map_rows.append(
            dict(
                id=r["id"],
                loadKg=r["loadKg"],
                speedCode=r["speedCode"],
                split=r["split"],
                rmse=float(np.sqrt(np.mean((r["mappedBack"] - r["back"]) ** 2))),
                correlation=float(np.corrcoef(r["mappedBack"], r["back"])[0, 1]),
                outOfDomain=bool(r["mapOutOfDomain"]),
            )
        )
    savemat(
        OUT / "map_model.mat",
        dict(
            mapModel=model,
            records=structs(
                [
                    {
                        k: v
                        for k, v in r.items()
                        if k not in ("profile", "six", "weights", "epoch")
                    }
                    for r in records
                ]
            ),
            trainingTrialPaths=np.array([r["trialPath"] for r in train], dtype=object),
        ),
        do_compression=True,
        long_field_names=True,
    )
    map_report = dict(
        subject=SUBJECT,
        architecture=model["architecture"],
        shapeFactor=0.1,
        shapeFactorChanged=False,
        nativeResponseAudits=shape_audits,
        trainingLoadsKg=sorted({r["loadKg"] for r in train}),
        trainingConditions=len(train),
        heldOutConditions=len(records) - len(train),
        metrics=map_rows,
        qualityReview=review,
        usesIKForMAP=False,
        usesTorqueForMAP=False,
    )
    (OUT / "map_evaluation.json").write_text(json.dumps(map_report, indent=2) + "\n")
    fig, axes = plt.subplots(5, 3, figsize=(13, 12), squeeze=False)
    for ax, r in zip(
        axes.ravel(), sorted(records, key=lambda r: (r["loadKg"], r["speedCode"]))
    ):
        ax.plot(
            np.arange(len(r["q"])) * 0.01, r["back"], label="Experimental activation"
        )
        ax.plot(np.arange(len(r["q"])) * 0.01, r["mappedBack"], label="Trained MAP")
        ax.set_title(
            f"{r['loadKg']:g} kg / code {r['speedCode']} | {r['split']}", fontsize=9
        )
        ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(
        "240617Pilot | excitation-only MAP training and whole-condition holdouts"
    )
    fig.tight_layout()
    fig.savefig(OUT / "map_training.png", dpi=150)
    plt.close(fig)

    choices = []
    audit = []
    for r in records:
        if r["loadKg"] not in (5, 10, 15):
            continue
        path = DATA / SUBJECT / r["trialPath"]
        try:
            tables, clocks, identity, exact, issues, repairs = native_tables(path)
            if issues:
                raise ValueError(str(issues))
            q, angle_audit = angle_reference(path, r["epoch"])
            if len(find_peaks(-q, prominence=0.12, distance=50)[0]) < 2:
                raise ValueError("Fewer than two measured flexion valleys")
            table = tables["Torque.sto"][1]
            valid = supported(table[:, 0], r["epoch"])
            if np.mean(valid) < 0.25:
                raise ValueError("Less than 25% torque reference support")
            tq = np.interp(r["epoch"], table[:, 0], table[:, 1])
            tq[~valid] = np.nan
            r = dict(
                r,
                q=q,
                torque=tq,
                torqueReferenceValid=valid,
                angleUsable=True,
                angleTorqueCorrelation=float(np.corrcoef(q[valid], tq[valid])[0, 1]),
                rightCensored=False,
                split="calibration",
                fitWeight=1.0,
                boundaryResidualWeight=2.0,
            )
            choices.append(r)
            audit.append(
                dict(
                    id=r["id"],
                    accepted=True,
                    referenceFraction=float(np.mean(valid)),
                    angle=angle_audit,
                    clocks=clocks,
                    torqueIdentity=identity,
                    appliedClockShiftSeconds=0.0,
                    rawSHA256={n: digest(path / n) for n in tables},
                )
            )
        except Exception as exc:
            audit.append(dict(id=r["id"], accepted=False, reason=str(exc)))
    # The speed-coded skeletal exports are often cut off mid-action. Retain
    # nine complete, previously audited bouts, with their actual cadence and
    # unknown nominal speed (-1). Repetitions are explicitly identified.
    quality_root = ROOT / "results/multi_subject_quality_v2"
    quality = json.loads((quality_root / "quality_gate.json").read_text())
    trials = {r["trial"]: r for r in quality["trials"] if r["subject"] == SUBJECT}
    selected = []
    qualified = sorted(many(old["cycles"]), key=lambda r: (r["loadKg"], r["id"]))
    for load, count in ((5, 5), (10, 2), (15, 2)):
        for original in [r for r in qualified if r["loadKg"] == load][:count]:
            trial = trials[original["trialPath"]]
            for name, expected in trial["fileSha256"].items():
                assert digest(DATA / SUBJECT / original["trialPath"] / name) == expected
            cache = quality_root / "signals" / trial["signalCache"]
            assert digest(cache) == trial["signalCacheSHA256"]
            signal = dict(np.load(cache))
            r = dict(original)
            clock = r["sourceStartSeconds"] + np.arange(len(r["q"])) * 0.01
            six = np.column_stack(
                [np.interp(clock, signal["t"], v) for v in signal["six"].T]
            )
            r["back"] = six @ model["lumpWeights"]
            r["mappedBack"] = map_predict(model, r)
            r["torqueReferenceValid"] = np.ones(len(r["q"]), dtype=bool)
            r["split"] = "calibration"
            r["fitWeight"] = 1.0
            r["boundaryResidualWeight"] = 2.0
            selected.append(r)
    assert len(selected) == 9
    generator = input_generator(
        [
            dict(
                a=np.c_[r["back"], r["mappedBack"]],
                front=r["front"],
                m=np.full(len(r["q"]), r["loadKg"]),
                loadKg=r["loadKg"],
                q=r["q"],
            )
            for r in train
        ]
    )
    groups = [
        dict(group=j + 1, sequence=combine([r], SUBJECT, r["id"]))
        for j, r in enumerate(selected)
    ]
    example = [next(r for r in selected if r["loadKg"] == load) for load in (5, 10, 15)]
    groups.append(
        dict(group=10, sequence=combine(example, SUBJECT, "Continuous 5/10/15 kg"))
    )
    old_models = loadmat(
        ROOT / "results/multi_subject_quality_v2" / SUBJECT / "neural_models.mat",
        simplify_cells=True,
    )
    models = [
        dict(m, generator=generator, intentDecoder="subtractive")
        for m in many(old_models["models"])
    ]
    savemat(
        OUT / "calibration_input.mat",
        dict(
            subject=SUBJECT,
            models=structs(models),
            controller=old_models["controller"],
            groups=structs(groups),
            mapModel=model,
        ),
        do_compression=True,
        long_field_names=True,
    )
    report = dict(
        subject=SUBJECT,
        selectionBeforeBSEFitting=True,
        conditions=[
            dict(
                id=r["id"],
                loadKg=r["loadKg"],
                speedCode=r["speedCode"],
                cadenceBPM=r["cadenceBPM"],
            )
            for r in selected
        ],
        nineDistinctNominalLoadSpeedPairs=False,
        composition="Five 5 kg bouts, two 10 kg bouts and two 15 kg bouts; measured cadence, unknown nominal speed",
        audit=audit,
        originalAngleScaleRetained=True,
        rawFilesChanged=False,
        usesIndependentTorque=False,
    )
    (OUT / "bse_selection.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                subject=SUBJECT,
                MAPconditions=len(records),
                MAPtrainingLoadsKg=map_report["trainingLoadsKg"],
                BSEconditions=report["conditions"],
            ),
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
