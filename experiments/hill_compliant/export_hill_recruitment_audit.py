"""Verify and export recruitment diagnostics without changing fitted models.

Run demo_hill_recruitment_audit first. The input mapping is reconstructed
independently here; successful numerical checks are not physiological validation.
"""
import csv
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.signal import lfilter

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/hill_recruitment_audit"
OUT = ROOT / "docs/hill_recruitment_audit"


def read_csv(name):
    with (RESULTS / name).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def memory(x, tau, dt):
    decay = np.exp(-dt / tau)
    return lfilter([1 - decay], [1, -decay], x, zi=[decay * x[0]])[0]


def residual(actual, expected, tolerance=1e-10):
    difference = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
    assert difference < tolerance, (difference, tolerance)
    return difference


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shade_waits(ax, trial):
    t = trial["t"]
    waiting = np.asarray(trial["waiting"], dtype=bool)
    edges = np.diff(np.r_[False, waiting, False].astype(int))
    for start, stop in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
        ax.axvspan(t[start], t[min(stop, len(t) - 1)], color=".94", zorder=0)
    ax.set_xlim(t[0], t[-1])
    ax.grid(alpha=.25)


def plot_summary(group, record):
    trial = record["baseline"]["sequence"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 7), sharex=True, sharey="row", layout="constrained")
    for source, title in enumerate(("Experimental activation", "MAP activation")):
        p = record["baseline"]["predictions"][source]
        base = record["details"][source]["base"]
        axes[0, source].plot(trial["t"], base[:, 1], color="#555555", label="Front base input")
        axes[0, source].plot(trial["t"], p["effectiveActivation"][:, 1], color="#dc6500", label="Abdominal Hill activation")
        axes[0, source].set_ylim(0, 1)
        axes[0, source].set_title(title, fontsize=12)
        axes[0, source].set_ylabel("Activation (model units; not verified %MVC)")
        axes[0, source].legend(loc="upper left", fontsize=9)
        axes[1, source].plot(trial["t"], p["backTorque"], color="#1676bb", label="Back extensor magnitude")
        axes[1, source].plot(trial["t"], p["frontTorque"], color="#dc6500", label="Abdominal flexor magnitude")
        axes[1, source].plot(trial["t"], p["torque"], "k--", label="Net Hill = back - abdomen")
        axes[1, source].set_ylabel("Torque (N·m)")
        axes[1, source].set_xlabel("Time (s)")
        axes[1, source].legend(loc="upper left", fontsize=9)
        for row in range(2):
            shade_waits(axes[row, source], trial)
    fig.suptitle(f"Group {group}: input evidence and opposing model torques | Gray: standing waits", fontsize=14)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"recruitment_summary_group{group}.{suffix}", dpi=180)
    plt.close(fig)


def main():
    audit = loadmat(RESULTS / "recruitment_audit.mat", simplify_cells=True)["audit"]
    sequence_path = ROOT / "results/hill_sequence/hill_sequence.mat"
    sequence = loadmat(sequence_path, simplify_cells=True)["result"]
    models_path = ROOT / "data/processed/hill_compliant/calibrated_models.mat"
    models = loadmat(models_path, simplify_cells=True)["models"]
    metrics = read_csv("recruitment_metrics.csv")
    variants = read_csv("counterfactual_metrics.csv")
    assert len(audit["groups"]) == len(sequence["groups"]) == 3
    assert len(metrics) == 6 and len(variants) == 24
    OUT.mkdir(parents=True, exist_ok=True)
    checks, components = [], []
    for group, record in enumerate(audit["groups"], 1):
        trial = record["baseline"]["sequence"]
        current = sequence["groups"][group - 1]
        residual(trial["q"], current["sequence"]["q"])
        residual(trial["torque"], current["sequence"]["torque"])
        moving = ~np.asarray(trial["waiting"], dtype=bool)
        dt = float(trial["dt"])
        for source in (0, 1):
            par = np.asarray(sequence["models"][source]["parameters"])
            residual(par, models[source]["parameters"])
            p = record["baseline"]["predictions"][source]
            d = record["details"][source]
            assert not p["failed"] and not np.any(d["failed"])
            for key in ("q", "v", "torque", "backTorque", "frontTorque", "effectiveActivation"):
                residual(p[key], current["predictions"][source][key])
            back = trial["a"][:, source]
            fast = memory(back, .1, dt)
            slow = memory(back, .35, dt)
            if source == 0:
                front = trial["front"]
            else:
                generator = sequence["models"][source]["generator"]
                m = trial["m"] / 15
                n = trial["loadKg"] / 15
                features = np.column_stack((back, back**2, fast, slow, back - slow, m, n,
                                            back * n, back * m, (back - fast) * n))
                features = (features - generator["mu"]) / generator["scale"]
                z = np.column_stack((np.ones(len(back)), features)) @ generator["beta"]
                front = 1 / (1 + np.exp(-np.clip(z, -30, 30)))
            base_residual = residual(d["base"], np.column_stack((back, front)))
            observed = par[14] * front
            inferred = par[15] * back + par[16] * np.maximum(0, fast - slow) + par[17] * np.maximum(0, slow - fast)
            recruitment = observed + inferred
            mapping_residual = residual(-np.expm1(-np.maximum(0, recruitment)), current["inputs"][source][:, 1])
            residual(d["observed"], observed)
            residual(d["inferred"], inferred)
            row = metrics[(group - 1) * 2 + source]
            calculated = {
                "frontBaseInputMean": np.mean(front[moving]),
                "frontHillActivationMean": np.mean(p["effectiveActivation"][moving, 1]),
                "frontHillActivationPeak": np.max(p["effectiveActivation"][moving, 1]),
                "backTorqueMeanNm": np.mean(p["backTorque"][moving]),
                "frontTorqueMeanNm": np.mean(p["frontTorque"][moving]),
                "frontToBackMeanTorqueRatio": np.sum(p["frontTorque"][moving]) / np.sum(p["backTorque"][moving]),
                "frontRecruitmentFromBackPercent": 100 * np.sum(inferred[moving]) / np.sum(recruitment[moving]),
            }
            for key, value in calculated.items():
                residual(float(row[key]), value)
            joint_residual = residual(p["torque"], p["backTorque"] - p["frontTorque"])
            component_residual = residual(p["frontTorque"], p["frontActive"] + p["frontPassive"] + p["frontFiberDamping"], 1e-7)
            trace = read_csv(f"group{group}_source{source + 1}.csv")
            for column, expected in (("abdominalTorqueNm", p["frontTorque"]),
                                     ("netHillTorqueNm", p["torque"]),
                                     ("frontBaseInput", front)):
                residual([float(r[column]) for r in trace], expected, 1e-9)
            for k, fraction in enumerate(d["fractions"]):
                variant = variants[(group - 1) * 8 + source * 4 + k]
                residual(float(variant["backDrivenFrontGainFraction"]), fraction)
                q = d["qVariants"][k]
                assert np.all(np.isfinite(q)) and float(variant["failed"]) == 0
                residual(float(variant["minimumAngleRad"]), np.min(q))
                residual(float(variant["maximumAngleRad"]), np.max(q))
                residual(float(variant["thetaRMSEDegrees"]), np.sqrt(np.mean((q[moving] - trial["q"][moving])**2)) * 180 / np.pi)
            peak = np.flatnonzero(moving)[np.argmax(p["frontTorque"][moving])]
            components.append({
                "group": group, "source": row["source"],
                "frontActiveMeanNm": float(np.mean(p["frontActive"][moving])),
                "frontPassiveMeanNm": float(np.mean(p["frontPassive"][moving])),
                "frontFiberDampingMeanNm": float(np.mean(p["frontFiberDamping"][moving])),
                "frontTorquePeakNm": float(p["frontTorque"][peak]),
                "frontActiveAtTotalPeakNm": float(p["frontActive"][peak]),
                "frontPassiveAtTotalPeakNm": float(p["frontPassive"][peak]),
            })
            checks.append({"group": group, "source": row["source"],
                           "baseInputResidual": base_residual, "inputMappingResidual": mapping_residual,
                           "netJointTorqueResidualNm": joint_residual,
                           "frontComponentResidualNm": component_residual})
        plot_summary(group, record)
    for name, records in (("numerical_audit.csv", checks), ("component_diagnostics.csv", components)):
        with (OUT / name).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    for name in ("recruitment_metrics.csv", "parameter_audit.csv", "counterfactual_metrics.csv"):
        shutil.copy2(RESULTS / name, OUT / name)
    for group in (1, 2, 3):
        for stem in ("recruitment_audit", "mapping_ablation"):
            for suffix in ("png", "pdf"):
                name = f"{stem}_group{group}.{suffix}"
                shutil.copy2(RESULTS / name, OUT / name)
    metadata = {
        "auditMatSHA256": digest(RESULTS / "recruitment_audit.mat"),
        "sequenceMatSHA256": digest(sequence_path), "modelsMatSHA256": digest(models_path),
        "physiologicalStatus": "Abdominal recruitment is not validated; large inferred cocontraction is unsupported by independent measurements.",
        "activationPolicy": "Project activation and Hill states are not established %MVC EMG.",
        "recruitmentPercentPolicy": "Fraction of pre-saturation mapping S from back/history, not a fraction of measured abdominal force or EMG.",
        "counterfactualPolicy": "Three input gains scaled with mechanics fixed; no parameter refit and no conclusion of physiological necessity.",
        "metricPolicy": "Movement frames only; means are sample-weighted, not a mean of per-condition means.",
        "files": [{"file": p.name, "sha256": digest(p)} for p in sorted(OUT.iterdir()) if p.suffix in (".csv", ".png", ".pdf")],
    }
    (OUT / "export_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("PASS: six input mappings, model/cache identity, opposing Nm torques and 24 counterfactual angle traces.")
    print("Numerical verification is not physiological validation.")
    for row in components:
        print(row)
    print("Snapshots:", OUT)


if __name__ == "__main__":
    main()
