"""Prepare subject-specific MAP calibration and continuous Hill replay data.

MAP retains four nonnegative synergies, three-Gaussian primitives and quadratic
load/cadence weight regression. Experimental angles keep the original project's
five-sample Gaussian smoothing and 0.5--1.4 scale; polarity is never flipped.
File/cycle holdouts are assigned before fitting. Missing conditions stay missing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat
from scipy.ndimage import binary_closing, gaussian_filter1d, label
from scipy.optimize import least_squares, nnls
from scipy.signal import find_peaks, lfilter

from audit_subjects import DATA, OUT, ROOT, extract

sys.path.insert(0, str(ROOT / "experiments/hill_compliant"))
from fit_hill_compliant import input_generator, pair_inputs


def structs(records):
    a = np.empty((1, len(records)), dtype=object)
    for i, r in enumerate(records):
        a[0, i] = r
    return a


def phase_features(load, cadence):
    return np.array([load / 15, (load / 15) ** 2, cadence / 10, (cadence / 10) ** 2])


def cycles(row, checked_signal=None):
    s = (
        checked_signal
        if checked_signal is not None
        else extract(DATA / row["subject"] / row["trial"])
    )
    dt = 0.01
    t = np.arange(0, s["t"][-1] + 1e-6, dt)
    signal = {
        key: np.interp(t, s["t"], s[key]) for key in ("q", "back", "front", "torque")
    }
    signal["six"] = np.column_stack([np.interp(t, s["t"], y) for y in s["six"].T])
    signal["weights"] = np.column_stack(
        [np.interp(t, s["t"], y) for y in s["weights"].T]
    )
    env = gaussian_filter1d(signal["back"], 12)
    baseline, peak = np.percentile(env, [10, 98])
    if peak - baseline < 0.07:
        return []
    active = env > baseline + 0.18 * (peak - baseline)
    original_active = active.copy()
    active = binary_closing(active, structure=np.ones(71))
    regions, n = label(active)
    result = []
    for region in range(1, n + 1):
        ids = np.flatnonzero(regions == region)
        if len(ids) < 90:
            continue
        left, right = max(0, ids[0] - 25), min(len(t) - 1, ids[-1] + 25)
        if right - left > 1600:
            continue
        # Start/end padding is from the recording, not an invented angle ramp.
        index = np.arange(left, right + 1)
        duration = float((right - left) * dt)
        q, tq = signal["q"][index], signal["torque"][index]
        qcorr = float(np.corrcoef(q, tq)[0, 1]) if np.std(q) > 1e-6 else 0.0
        valleys = find_peaks(-q, prominence=0.12, distance=50)[0]
        profile = np.column_stack(
            [
                np.interp(
                    np.linspace(0, duration, 200),
                    t[index] - t[left],
                    signal["six"][index, j],
                )
                for j in range(6)
            ]
        )
        result.append(
            {
                "id": row["trial"].replace("\\", "_") + f"_cycle{region}",
                "trialPath": row["trial"],
                "loadKg": row["loadKg"],
                "speedCode": row["speedCode"],
                "cadenceBPM": 60 / duration,
                "duration": duration,
                "sourceStartSeconds": float(t[left]),
                "sourceEndSeconds": float(t[right]),
                "q": q,
                "torque": tq,
                "back": signal["back"][index],
                "front": signal["front"][index],
                "six": signal["six"][index],
                "weights": signal["weights"][index],
                "profile": profile,
                "angleTorqueCorrelation": qcorr,
                "angleUsable": bool(
                    qcorr < 0.15
                    and np.ptp(q) > 0.2
                    and (
                        checked_signal is not None
                        or (
                            row["subject"] != "240723Pilot"
                            and not (
                                row["subject"] == "240719Pilot" and row["loadKg"] == 10
                            )
                        )
                    )
                ),
                "referenceValleys": valleys,
                "split": "train",
            }
        )
        result[-1]["rightCensored"] = bool(original_active[-1] and right >= len(t) - 70)
    return result


def assign_splits(subject, records):
    paths = sorted(set(r["trialPath"] for r in records))
    for path in paths:
        group = [r for r in records if r["trialPath"] == path]
        name = Path(path).name
        whole_test = (
            (subject == "240617Pilot" and name.endswith("_33"))
            or (
                subject == "240627Pilot"
                and (
                    name.endswith("_45_1")
                    and "5kg" not in name
                    or name.endswith("_45_2")
                )
            )
            or (subject == "240723Pilot" and name.endswith("_old"))
            or (subject == "240725Pilot" and name.endswith("5kg_3"))
        )
        if whole_test:
            for r in group:
                r["split"] = "test_file"
        elif len(group) >= 3:
            group[-1]["split"] = "test_cycle"
    return records


def qualified_splits(subject, records):
    """Keep a holdout per load even when the old holdout file fails quality."""
    records = assign_splits(subject, records)
    # In the historical splitter '5kg' was also a substring of '15kg'. Keep
    # historical frozen replay stable, and correct this only for new fits.
    if subject == "240627Pilot":
        for r in records:
            if r["loadKg"] == 15 and Path(r["trialPath"]).name.endswith("_45_1"):
                r["split"] = "test_file"
    for load in sorted(set(r["loadKg"] for r in records)):
        group = [r for r in records if r["loadKg"] == load]
        if len(group) >= 2 and not any(r["split"] != "train" for r in group):
            group[-1]["split"] = "test_qualified_cycle"
        if group and not any(r["split"] == "train" for r in group):
            # Preserve at least one training record; a single qualified action
            # is calibration-only and cannot validate that load independently.
            group[0]["split"] = "train"
    return records


def gaussian_model(z, phase):
    z = z.reshape(3, 3)
    return sum(a * np.exp(-(((phase - c) / w) ** 2)) for a, c, w in z)


def fit_map(train):
    y = np.hstack([r["profile"] for r in train])
    rng = np.random.default_rng(617)
    w, h = rng.uniform(0.1, 0.8, (200, 4)), rng.uniform(0.01, 0.3, (4, y.shape[1]))
    for _ in range(1200):
        h *= (w.T @ y) / np.maximum(w.T @ w @ h + 0.0001 * h, 1e-12)
        w *= (y @ h.T) / np.maximum(w @ h @ h.T + 0.0001 * w, 1e-12)
        scale = np.maximum(np.max(w, axis=0), 1e-10)
        w /= scale
        h *= scale[:, None]
    phase = np.linspace(0, 1, 200)
    gaussian_parameters = []
    for j in range(4):
        best = None
        for centers in ((0.2, 0.5, 0.8), (0.3, 0.45, 0.7), (0.15, 0.65, 0.85)):
            z = np.array([[0.6, c, 0.15] for c in centers]).ravel()
            fit = least_squares(
                lambda z: gaussian_model(z, phase) - w[:, j],
                z,
                bounds=(np.tile([-2.0, 0.0, 0.025], 3), np.tile([3.0, 1.0, 0.6], 3)),
                max_nfev=500,
            )
            if best is None or fit.cost < best.cost:
                best = fit
        gaussian_parameters.append(best.x)
    primitives = np.maximum(
        0, np.column_stack([gaussian_model(z, phase) for z in gaussian_parameters])
    )
    target = np.array(
        [
            np.column_stack([nnls(primitives, r["profile"][:, j])[0] for j in range(6)])
            for r in train
        ]
    )
    x = np.array([phase_features(r["loadKg"], r["cadenceBPM"]) for r in train])
    mu = x.mean(0)
    scale = np.maximum(x.std(0), 0.1)
    design = np.c_[np.ones(len(x)), (x - mu) / scale]
    penalty = np.eye(5) * 0.15
    penalty[0, 0] = 1e-8
    beta = np.linalg.solve(
        design.T @ design + penalty, design.T @ target.reshape(len(train), -1)
    )
    # Fixed training weights prevent MAP from consuming future geometry/torque.
    lump_weights = np.mean(np.vstack([r["weights"] for r in train]), axis=0)
    lump_weights /= lump_weights.sum()
    return {
        "primitives": primitives,
        "gaussianParameters": np.array(gaussian_parameters),
        "mu": mu,
        "scale": scale,
        "covariateMinimum": x.min(0),
        "covariateMaximum": x.max(0),
        "beta": beta,
        "lumpWeights": lump_weights,
        "phase": phase,
        "architecture": "4 NMF synergies; 3-Gaussian primitives; quadratic load/cadence weight regression",
    }


def map_predict(model, record):
    features = phase_features(record["loadKg"], record["cadenceBPM"])
    bounded = np.clip(features, model["covariateMinimum"], model["covariateMaximum"])
    record["mapOutOfDomain"] = bool(np.any(np.abs(features - bounded) > 1e-10))
    design = np.r_[1.0, (bounded - model["mu"]) / model["scale"]]
    weights = np.maximum(0, design @ model["beta"]).reshape(4, 6)
    profile = model["primitives"] @ weights
    x = np.linspace(0, 1, len(record["q"]))
    six = np.column_stack(
        [np.interp(x, model["phase"], profile[:, j]) for j in range(6)]
    )
    return np.clip(six @ model["lumpWeights"], 0, 1)


def combine(records, subject, name):
    # One continuously evolving state through all recorded actions and waits.
    arrays = {
        k: []
        for k in (
            "q",
            "torque",
            "a",
            "front",
            "m",
            "loadKg",
            "speedCode",
            "stageId",
            "waiting",
            "referenceIsFilled",
        )
    }
    partial_torque = any("torqueReferenceValid" in r for r in records)
    if partial_torque:
        arrays["torqueReferenceValid"] = []
    segments = []
    offset = 0

    def wait(seconds):
        nonlocal offset
        n = round(seconds / 0.01)
        for key in arrays:
            value = np.full(
                n,
                (
                    1.5
                    if key == "q"
                    else (
                        20.0
                        if key == "torque"
                        else 1.0 if key in ("waiting", "referenceIsFilled") else 0.0
                    )
                ),
            )
            if key == "a":
                value = np.zeros((n, 2))
            arrays[key].append(value)
        offset += n

    wait(2.0)
    for stage, r in enumerate(records, 1):
        n = len(r["q"])
        start = offset
        a = np.c_[r["back"], r["mappedBack"]]
        mass = np.full(n, r["loadKg"])
        # Contact is a declared schedule inferred from the EMG cycle boundaries,
        # not a measured force sensor. Release is smoothed for 80 ms.
        raw_mass = np.r_[
            np.zeros(min(25, n)),
            np.full(max(0, n - 50), r["loadKg"]),
            np.zeros(min(25, n)),
        ][:n]
        if "loadContactStartSeconds" in r:
            clock=np.arange(n)*.01
            raw_mass=r["loadKg"]*((clock>=r["loadContactStartSeconds"]) &
                                   (clock<r["loadContactEndSeconds"])).astype(float)
        decay = np.exp(-0.01 / 0.08)
        mass = lfilter(1 - decay, [1, -decay], raw_mass)
        values = {
            "q": r["q"],
            "torque": r["torque"],
            "a": a,
            "front": r["front"],
            "m": mass,
            "loadKg": np.full(n, r["loadKg"]),
            "speedCode": np.full(n, r["speedCode"]),
            "stageId": np.full(n, stage),
            "waiting": np.zeros(n),
            "referenceIsFilled": np.zeros(n),
        }
        if partial_torque:
            values["torqueReferenceValid"] = np.asarray(
                r.get("torqueReferenceValid", np.isfinite(r["torque"])), dtype=bool
            )
        for key, value in values.items():
            arrays[key].append(value)
        offset += n
        segments.append(
            {
                "loadKg": r["loadKg"],
                "speedCode": r["speedCode"],
                "cadenceBPM": r["cadenceBPM"],
                "startIndex": start + 1,
                "endIndex": offset,
                "startTime": start * 0.01,
                "endTime": (offset - 1) * 0.01,
                "recordId": r["id"],
                "split": r["split"],
                "angleUsable": r["angleUsable"],
                "angleTorqueCorrelation": r["angleTorqueCorrelation"],
                **{k: r[k] for k in ("fitWeight", "boundaryResidualWeight", "loadContactStartSeconds", "loadContactEndSeconds") if k in r},
            }
        )
        wait(2.5)
    data = {key: np.concatenate(value) for key, value in arrays.items()}
    decay = np.exp(-0.01 / 0.12)
    data["front"] = lfilter(1 - decay, [1, -decay], data["front"])
    return {
        **data,
        "name": name,
        "subject": subject,
        "t": np.arange(offset) * 0.01,
        "dt": 0.01,
        "standingAngle": 1.5,
        "v0": 0.0,
        "segments": structs(segments),
        "referencePolicy": "Original Gaussian-5 and min-max 0.5--1.4 per recording; filled waits 1.5 rad/20 Nm explicitly flagged",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen-map-directory", type=Path)
    parser.add_argument("--allow-training-changes", action="store_true")
    parser.add_argument("--quality-gate", type=Path)
    parser.add_argument("--subjects", nargs="+")
    parser.add_argument("--loads", nargs="+", type=float, default=[5, 10, 15])
    parser.add_argument("--held-out-trials", nargs="+", default=[])
    parser.add_argument("--map-model-directory", type=Path)
    args = parser.parse_args()
    if args.allow_training_changes and args.quality_gate is None:
        parser.error("Recalibration requires --quality-gate from quality_gate.py.")
    gate = None
    if args.quality_gate:
        gate = json.loads(args.quality_gate.read_text(encoding="utf-8"))
        if gate["version"] != "clock-smoothing-gate-v2":
            raise RuntimeError("Unsupported quality gate version")
    audit = json.loads((OUT / "data_quality.json").read_text())
    for subject in sorted(set(r["subject"] for r in audit)):
        if args.subjects and subject not in args.subjects:
            continue
        rows = [
            r
            for r in audit
            if r["subject"] == subject
            and r["status"] in ("usable", "review")
            and r["loadKg"] in args.loads
            and r.get("backMax", 0) > 0.15
        ]
        records = []
        for row in rows:
            if gate is None:
                records.extend(cycles(row))
                continue
            checked = next(
                (
                    r
                    for r in gate["trials"]
                    if r["subject"] == subject and r["trial"] == row["trial"]
                ),
                None,
            )
            if checked is None or checked["gate"] != "pass":
                continue
            for name, expected in checked["fileSha256"].items():
                raw = DATA / subject / row["trial"] / name
                if hashlib.sha256(raw.read_bytes()).hexdigest() != expected:
                    raise RuntimeError(
                        "Raw files changed after quality audit: " + str(raw)
                    )
            signal_path = args.quality_gate.parent / "signals" / checked["signalCache"]
            if (
                hashlib.sha256(signal_path.read_bytes()).hexdigest()
                != checked["signalCacheSHA256"]
            ):
                raise RuntimeError("Audited signal cache changed: " + str(signal_path))
            signal = dict(np.load(signal_path))
            accepted = {
                r["id"]
                for r in gate["cycles"]
                if r["subject"] == subject and r["accepted"]
            }
            records.extend(r for r in cycles(row, signal) if r["id"] in accepted)
        records = (
            qualified_splits(subject, records)
            if gate is not None
            else assign_splits(subject, records)
        )
        for r in records:
            if Path(r["trialPath"]).name in args.held_out_trials:
                r["split"] = "test_file"
        if not records:
            continue
        train = [r for r in records if r["split"] == "train"]
        if gate is not None and (
            len(train) < 3
            or len(records) - len(train) < 2
            or len(set(r["loadKg"] for r in train)) < 2
        ):
            print(
                subject,
                "NOT CALIBRATED: insufficient qualified training/holdout/load coverage",
                flush=True,
            )
            continue
        frozen = None
        if args.frozen_map_directory is not None:
            frozen = loadmat(
                args.frozen_map_directory / (subject + "_models.mat"),
                simplify_cells=True,
            )
            model = frozen["mapModel"]
        elif args.map_model_directory is not None:
            external = loadmat(
                args.map_model_directory / (subject + "_map.mat"), simplify_cells=True
            )
            held_paths = {r["trialPath"] for r in records if r["split"] != "train"}
            map_paths = set(np.atleast_1d(external["trainingTrialPaths"]).tolist())
            if held_paths & map_paths:
                raise RuntimeError("MAP training overlaps BSE held-out trials")
            model = external["mapModel"]
        else:
            model = fit_map(train)
        metrics = []
        for r in records:
            r["mappedBack"] = map_predict(model, r)
            rmse = np.sqrt(np.mean((r["mappedBack"] - r["back"]) ** 2))
            metrics.append(
                {
                    k: r[k]
                    for k in (
                        "id",
                        "trialPath",
                        "loadKg",
                        "speedCode",
                        "cadenceBPM",
                        "duration",
                        "split",
                        "angleUsable",
                        "angleTorqueCorrelation",
                        "rightCensored",
                        "mapOutOfDomain",
                    )
                }
            )
            metrics[-1].update(
                mapRMSE=float(rmse),
                mapCorrelation=float(np.corrcoef(r["mappedBack"], r["back"])[0, 1]),
            )
        # Fit the existing causal abdominal generator on training activation only.
        gen_trials = [
            {
                "a": np.c_[r["back"], r["mappedBack"]],
                "front": r["front"],
                "m": np.full(len(r["q"]), r["loadKg"]),
                "loadKg": r["loadKg"],
                "q": r["q"],
            }
            for r in train
        ]
        generator = (
            np.asarray(frozen["models"], dtype=object).reshape(-1)[0]["generator"]
            if frozen is not None
            else input_generator(gen_trials)
        )
        # Separate sequences per recorded trial keep holdout provenance visible.
        trial_sequences = [
            combine([r for r in records if r["trialPath"] == p], subject, p)
            for p in sorted(set(r["trialPath"] for r in records))
        ]
        train_sequences = [combine([r], subject, r["id"]) for r in train]
        path = OUT / subject
        path.mkdir(parents=True, exist_ok=True)
        if not args.allow_training_changes and (path / "prepared.mat").exists():
            old = loadmat(path / "prepared.mat", simplify_cells=True)
            previous = np.asarray(old["trainSequences"], dtype=object).reshape(-1)
            assert len(previous) == len(train_sequences)
            for first, second in zip(previous, train_sequences):
                assert np.allclose(
                    first["a"], second["a"], rtol=0, atol=1e-12
                ), "Changing training inputs requires model recalibration."
        savemat(
            path / "prepared.mat",
            {
                "sequences": structs(trial_sequences),
                "trainSequences": structs(train_sequences),
                "cycles": structs(
                    [
                        {
                            k: v
                            for k, v in r.items()
                            if k not in ("profile", "six", "weights")
                        }
                        for r in records
                    ]
                ),
                "mapModel": model,
                "generator": generator,
                "qualityGateVersion": (
                    gate["version"] if gate else "historical_unchecked_replay"
                ),
                "qualityGateSHA256": (
                    hashlib.sha256(args.quality_gate.read_bytes()).hexdigest()
                    if gate
                    else ""
                ),
            },
            do_compression=True,
            long_field_names=True,
        )
        (path / "preparation.json").write_text(
            json.dumps(
                {
                    "subject": subject,
                    "cycles": metrics,
                    "trainingCycles": len(train),
                    "testCycles": len(records) - len(train),
                    "angleUsableCycles": sum(r["angleUsable"] for r in records),
                    "knownConditions": sorted(
                        set((r["loadKg"], r["speedCode"]) for r in records)
                    ),
                    "mapArchitecture": model["architecture"],
                    "mapUsesReferenceAngle": False,
                    "covariatesBoundedToTrainingDomain": True,
                    "loadContactIsMeasured": False,
                    "unknownSpeedCode": -1,
                    "qualityGateVersion": (
                        gate["version"] if gate else "historical_unchecked_replay"
                    ),
                    "qualityGateSHA256": (
                        hashlib.sha256(args.quality_gate.read_bytes()).hexdigest()
                        if gate
                        else ""
                    ),
                    "qualityCheckedBeforeTraining": gate is not None,
                    "explicitHeldOutTrials": args.held_out_trials,
                    "externalMapModelSHA256": (
                        hashlib.sha256(
                            (
                                args.map_model_directory / (subject + "_map.mat")
                            ).read_bytes()
                        ).hexdigest()
                        if args.map_model_directory
                        else ""
                    ),
                    "activationSmoothing": (
                        gate["filterPolicy"] if gate else "historical"
                    ),
                },
                indent=2,
            )
            + "\n"
        )
        print(
            subject,
            "cycles",
            len(records),
            "train",
            len(train),
            "angle usable",
            sum(r["angleUsable"] for r in records),
            "MAP train/test RMSE",
            [
                (
                    split,
                    float(
                        np.mean(
                            [
                                m["mapRMSE"]
                                for m in metrics
                                if (m["split"] == "train") == split
                            ]
                        )
                    ),
                )
                for split in (True, False)
            ],
            flush=True,
        )


if __name__ == "__main__":
    main()
