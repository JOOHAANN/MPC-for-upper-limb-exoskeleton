"""Reduce unnecessary support while keeping back torque just below its target.

Only qualified training actions select controller options. BSE/MAP parameters,
references and activation streams stay byte-identical to the quality-first run.
The unassisted trajectory supplies a demand gate, not a measured reference.
"""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.io import loadmat
from scipy.signal import find_peaks

from audit_subjects import ROOT
from prepare_subjects import combine
from fit_subject import many, record, rollout
from control_subject import metrics
from run_hill_mpc import simulate_case

BASE = ROOT / "results/multi_subject_quality_v2"
DEST = ROOT / "results/mpc_threshold_tracking"


def setup(subject):
    output = DEST / subject
    output.mkdir(parents=True, exist_ok=True)
    for name in (
        "prepared.mat",
        "neural_models.mat",
        "hill_sequence.mat",
        "preparation.json",
        "fit_policy.json",
        "bse_metrics.json",
    ):
        source = BASE / subject / name
        target = output / name
        if not target.exists() or target.read_bytes() != source.read_bytes():
            shutil.copy2(source, target)
    for name in ("quality_gate.json", "data_quality.json"):
        if (BASE / name).exists():
            shutil.copy2(BASE / name, DEST / name)
    return output


def tracking_metrics(actual, nominal, seq, source):
    rows = metrics(actual, nominal, seq, np.zeros(len(actual)), source)
    for row, seg in zip(rows, many(seq["segments"])):
        ids = np.arange(int(seg["startIndex"]) - 1, int(seg["endIndex"]))
        limit = row["backThresholdNm"]
        demand = nominal[ids, 3] > limit
        row.update(
            demandSeconds=float(sum(demand) * seq["dt"]),
            demandMeanShortfallNm=(
                float(np.mean(np.maximum(0, limit - actual[ids[demand], 3])))
                if np.any(demand)
                else 0.0
            ),
            demandMeanAssistanceNm=(
                float(np.mean(actual[ids[demand], 24])) if np.any(demand) else 0.0
            ),
            assistanceImpulseNms=float(np.sum(actual[ids, 24]) * seq["dt"]),
            nominalPeakExceedsTarget=bool(nominal[ids, 3].max() > limit),
            peakShortfallNm=(
                float(max(0, limit - actual[ids, 3].max()))
                if nominal[ids, 3].max() > limit
                else None
            ),
        )
    return rows


def score(rows):
    demand = sum(r["demandSeconds"] for r in rows)
    excess = max(r["excessPeakNm"] for r in rows)
    valleys = sum(abs(r["predictedValleys"] - r["nominalValleys"]) for r in rows)
    deficit = sum(r["demandMeanShortfallNm"] * r["demandSeconds"] for r in rows) / max(
        demand, 1e-9
    )
    impulse = sum(r["assistanceImpulseNms"] for r in rows)
    return (excess > 0.1, excess if excess > 0.1 else 0.0, valleys, deficit, impulse)


def main():
    global BASE, DEST
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("subject")
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--destination", type=Path, default=DEST)
    parser.add_argument("--scales", nargs="+", type=float, default=[0.25, 0.5, 0.75])
    parser.add_argument("--under-weight", type=float, default=5.0)
    parser.add_argument("--effort-weight", type=float, default=4.0)
    parser.add_argument("--guard-floor", type=float, default=3.0)
    parser.add_argument("--horizon", type=float)
    parser.add_argument(
        "--screen",
        action="store_true",
        help="Screen one training action per load before a full training pass.",
    )
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    BASE, DEST = args.base, args.destination
    if any(not 0.0 <= value <= 1.0 for value in args.scales):
        parser.error("Guard scales must be between zero and one.")
    if min(args.under_weight, args.effort_weight, args.guard_floor) < 0:
        parser.error("Weights and guard floor must be nonnegative.")
    if args.horizon is not None and not 0.025 <= args.horizon <= 0.4:
        parser.error("The constraint horizon must be between 0.025 and 0.4 seconds.")
    output = setup(args.subject)
    data = loadmat(output / "prepared.mat", simplify_cells=True)
    assert data["qualityGateVersion"] in (
        "clock-smoothing-gate-v2",
        "excitation-ik-three-loads-v1",
    )
    cache = loadmat(output / "neural_models.mat", simplify_cells=True)
    models = many(cache["models"])
    cfg = cache["controller"]
    original_file = BASE / args.subject / "mpc_evaluation.json"
    if not original_file.exists():
        original_file = BASE / args.subject / "mpc_selected_training_options.json"
    original = json.loads(original_file.read_text())["options"]
    train = [c for c in many(data["cycles"]) if c["split"] == "train"]
    if args.screen:
        train = [
            next(c for c in train if c["loadKg"] == load)
            for load in sorted(set(c["loadKg"] for c in train))
        ]
    seq = combine(train, args.subject, "threshold_tracking_training")
    rec = record(seq, data["generator"])
    inputs = [
        rollout(m["parameters"], m["neuralParameters"], cfg, rec, s, 10)[2]
        for s, m in enumerate(models)
    ]
    name = "screen_candidates.json" if args.screen else "threshold_candidates.json"
    audited = (
        json.loads((output / name).read_text())
        if args.resume and (output / name).exists()
        else []
    )
    notified = [time.monotonic()]

    def progress(k, n):
        if time.monotonic() - notified[0] > 25:
            print(args.subject, "TRACKING", round(k / n * 100, 1), "%", flush=True)
            notified[0] = time.monotonic()

    for scale in args.scales:
        opt = {
            **original,
            "predictionErrorGuardScale": scale,
            "predictionErrorGuardFloorNm": args.guard_floor,
            "backTorqueUnderWeight": args.under_weight,
            "trackingHeadroomNm": 3.0,
            "assistanceEffortWeight": args.effort_weight,
        }
        if args.horizon is not None:
            opt["torqueConstraintHorizonSeconds"] = args.horizon
        if any(c["options"] == opt for c in audited):
            continue
        rows = []
        for source in (0, 1):
            a, n, f, _ = simulate_case(
                models[source], cfg, seq, inputs[source], opt, progress
            )
            rows.extend(tracking_metrics(a, n, seq, source))
            assert a[:, 24].min() >= -1e-8 and a[:, 24].max() <= 80 + 1e-8
            assert max(abs(np.diff(a[:, 24]))) <= 2 + 1e-8
        audit = dict(
            options=opt,
            cases=rows,
            selectionScore=score(rows),
            trainingRecordIds=[r["id"] for r in train],
            modelSHA256=hashlib.sha256(
                (output / "neural_models.mat").read_bytes()
            ).hexdigest(),
            trainingOnly=True,
            screeningOnly=bool(args.screen),
        )
        audited.append(audit)
        (output / name).write_text(json.dumps(audited, indent=2) + "\n")
        print(
            args.subject,
            "CANDIDATE scale",
            scale,
            "SCORE",
            score(rows),
            "peaks",
            [round(r["backPeakNm"], 2) for r in rows],
            flush=True,
        )
    chosen = min(audited, key=lambda c: score(c["cases"]))
    if not args.screen:
        (output / "selected_threshold_options.json").write_text(
            json.dumps(chosen, indent=2) + "\n"
        )
    print(
        "SELECTED",
        args.subject,
        chosen["options"],
        "SCORE",
        score(chosen["cases"]),
        flush=True,
    )


if __name__ == "__main__":
    main()
