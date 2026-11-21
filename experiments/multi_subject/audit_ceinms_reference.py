"""Audit reference torque and export A=0.1 model candidates on native geometry.

No independent force/ID target exists for these Exp records. Candidates are
explicitly uncalibrated. State reconstruction is a same-state sensitivity
calculation, not an independent prediction or a replay of the unavailable DLL.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat

from audit_subjects import DATA, ROOT, read_sto
from ceinms_force_replay import StiffMuscles, causal_velocity
from compare_exp_signal_levels import activation_curve
from recompute_exp_activation import write_sto
from review_exp_one_cycle import sha

OUTPUT = ROOT / "docs/ceinms_reference_audit"
DERIVED = ROOT / "data/processed/ceinms_reference_audit"
REVIEW = ROOT / "docs/exp_one_cycle_review/review.json"
TABLES = (
    "MuscleForces",
    "MuscleActiveForces",
    "MusclePassiveForces",
    "FibreLengths",
    "FibreVelocities",
    "PennationAngle",
    "lmt",
    "ma_L5_S1_Flex_Ext",
    "Torque",
)


def replay(row, xml, state_export_directory=None):
    path = DATA / row["subject"] / row["trial"]
    tables = {name: read_sto(path / (name + ".sto")) for name in TABLES}
    common = tables["MuscleForces"][1][:, 0]
    for _, values in tables.values():
        common = np.intersect1d(common, values[:, 0])
    if len(common) < 30:
        raise ValueError("Insufficient common muscle-state samples")
    names = tables["MuscleForces"][0][1:]
    values = {}
    for name, (columns, table) in tables.items():
        indices = np.searchsorted(table[:, 0], common)
        values[name] = (
            table[indices, 1:]
            if name == "Torque"
            else table[indices][:, [columns.index(n) for n in names]]
        )
    model = StiffMuscles(xml, names)
    activation_names, activation = read_sto(path / "Activations.sto")
    initial = np.column_stack(
        [
            np.interp(
                common, activation[:, 0], activation[:, activation_names.index(n)]
            )
            for n in names
        ]
    )
    a0 = model.recover_recorded_activation(
        values["lmt"], values["FibreLengths"], values["MuscleActiveForces"], initial
    )
    old = model.forces(a0, values["FibreLengths"], np.zeros_like(a0))
    arms = values["ma_L5_S1_Flex_Ext"]
    torque = lambda force: np.sum(force * arms, axis=1)
    rmse = lambda data: float(np.sqrt(np.mean(data**2)))
    checks = dict(
        totalEqualsActiveOnlyRMSE_N=rmse(
            values["MuscleForces"] - np.maximum(0, values["MuscleActiveForces"])
        ),
        elasticMatchesPassiveExportRMSE_N=rmse(
            old["elastic"] - values["MusclePassiveForces"]
        ),
        reconstructedStaticActiveRMSE_N=rmse(
            old["active"] - values["MuscleActiveForces"]
        ),
        rigidLengthReplayRMSE_m=rmse(
            model.lengths(values["lmt"], a0) - values["FibreLengths"]
        ),
        forceArmTorqueIdentityRMSE_Nm=rmse(
            torque(values["MuscleForces"]) - values["Torque"][:, 0]
        ),
        causalVelocityReplayRMSE_m_per_s=rmse(
            causal_velocity(common, values["FibreLengths"]) - values["FibreVelocities"]
        ),
        standardFullFormulaVersusExportRMSE_N=rmse(
            model.forces(a0, values["FibreLengths"], values["FibreVelocities"])["total"]
            - values["MuscleForces"]
        ),
    )
    verified = (
        checks["totalEqualsActiveOnlyRMSE_N"] < 1e-5
        and checks["elasticMatchesPassiveExportRMSE_N"] < 1e-5
        and checks["reconstructedStaticActiveRMSE_N"] < 0.001
        and checks["rigidLengthReplayRMSE_m"] < 1e-6
    )
    record = dict(
        subject=row["subject"],
        trial=row["trial"],
        samples=len(common),
        loadKg=row["loadKg"],
        speedCode=row["speedCode"],
        checks=checks,
        originalEquivalentPolicyVerified=verified,
        originalConfiguration=str(xml),
        originalConfigurationSHA256=sha(xml),
        sourceSHA256={
            name + ".sto": sha(path / (name + ".sto"))
            for name in (*TABLES, "Activations")
        },
    )
    if not verified:
        record["candidateExported"] = False
        return record, None
    # Recover u from the OLD mapping, then evaluate the requested NEW mapping.
    # This avoids fitting a time shift between two asynchronous export streams.
    logarithm = 1 + a0 * np.expm1(model.old_shape)
    if np.min(logarithm) <= 0:
        raise ValueError("Cannot invert old activation mapping on recorded states")
    u = np.log(logarithm) / model.old_shape
    a1 = activation_curve(u, 0.1)
    lm1 = model.lengths(values["lmt"], a1)
    v1 = values["FibreVelocities"] + causal_velocity(
        common, lm1 - values["FibreLengths"]
    )
    static = model.forces(a1, lm1, np.zeros_like(v1))
    full = model.forces(a1, lm1, v1)
    back = np.array([n.startswith(("IL_", "LTpT_", "LTpL_", "MF_")) for n in names])
    abdomen = np.array([n.startswith(("rect_abd_", "EO", "IO")) for n in names])
    component = lambda force, mask: np.sum(force[:, mask] * arms[:, mask], axis=1)
    signals = dict(
        timeEpochSeconds=common,
        timeSeconds=common - row["originEpochSeconds"],
        originalTorque=values["Torque"][:, 0],
        A01OriginalPolicy=torque(static["activeOnly"]),
        A01PlusElastic=torque(static["total"]),
        A01FullFormulaCandidate=torque(full["total"]),
        oldPassiveBackTorque=component(values["MusclePassiveForces"], back),
        oldPassiveAbdomenTorque=component(values["MusclePassiveForces"], abdomen),
        oldPassiveOtherTorque=component(
            values["MusclePassiveForces"], ~(back | abdomen)
        ),
        A01FullBackTorque=component(full["total"], back),
        A01FullAbdomenTorque=component(full["total"], abdomen),
        A01FullOtherTorque=component(full["total"], ~(back | abdomen)),
        shapeFactor=0.1,
        externallyCalibrated=False,
    )
    destination = DERIVED / row["subject"] / row["trial"]
    destination.mkdir(parents=True, exist_ok=True)
    for key, filename in (
        ("A01OriginalPolicy", "Torque_A01_original_policy.sto"),
        ("A01PlusElastic", "Torque_A01_plus_elastic_candidate.sto"),
        ("A01FullFormulaCandidate", "Torque_A01_full_formula_candidate.sto"),
    ):
        write_sto(
            destination / filename,
            ["time", "L5_S1_Flex_Ext"],
            np.column_stack([common, signals[key]]),
            "UncalibratedModelEstimate_" + key,
        )
    savemat(destination / "torque_sensitivity.mat", signals, do_compression=True)
    if state_export_directory is not None:
        state_export_directory.mkdir(parents=True, exist_ok=True)
        sine = model.l0 * np.sin(np.deg2rad(model.pennation_parameter)) / lm1
        pennation = np.rad2deg(np.arcsin(np.clip(sine, 0, 1)))
        states = {
            "MuscleForces": static["activeOnly"],
            "MuscleActiveForces": static["active"],
            "MusclePassiveForces": static["elastic"],
            "FibreLengths": lm1,
            "FibreVelocities": v1,
            "PennationAngle": pennation,
            "TendonLength": np.broadcast_to(model.slack, lm1.shape),
            "TendonStrain": np.zeros_like(lm1),
            "skeletal_activation_A01": a1,
            "skeletal_excitation_recovered": u,
        }
        for name, state in states.items():
            write_sto(
                state_export_directory / (name + ".sto"),
                ["time", *names],
                np.column_stack([common, state]),
                "A01_original_equivalent_policy_" + name,
            )
        write_sto(
            state_export_directory / "Torque.sto",
            ["time", "L5_S1_Flex_Ext"],
            np.column_stack([common, signals["A01OriginalPolicy"]]),
            "ModelEstimate_A01_original_equivalent_active_only_NOT_measured_truth",
        )
    record.update(
        candidateExported=True,
        referenceReplacementAuthorizedByEvidence=False,
        newActivationOutsidePhysicalRangeFraction=float(
            np.mean((a1 < -1e-6) | (a1 > 1 + 1e-6))
        ),
        passiveTorqueMinimumNm=float(torque(values["MusclePassiveForces"]).min()),
        passiveTorqueMaximumNm=float(torque(values["MusclePassiveForces"]).max()),
        oldPeakNm=float(signals["originalTorque"].max()),
        A01OriginalPolicyPeakNm=float(signals["A01OriginalPolicy"].max()),
        A01FullFormulaPeakNm=float(signals["A01FullFormulaCandidate"].max()),
        A01FullFormulaMinimumNm=float(signals["A01FullFormulaCandidate"].min()),
        derivedDirectory=str(destination),
    )
    return record, signals


def examples(review, records, signals):
    keys = (("240725Exp", 5, 45), ("240730Exp", 10, 45))
    fig, axes = plt.subplots(3, 2, figsize=(12, 9))
    data = []
    for j, key in enumerate(keys):
        row = next(
            r
            for r in review["conditions"]
            if (r["subject"], r["loadKg"], r["speedCode"]) == key
        )
        s = signals[(row["subject"], row["trial"])]
        choice = row["selected"]
        keep = (s["timeEpochSeconds"] >= choice["startEpochSeconds"]) & (
            s["timeEpochSeconds"] <= choice["endEpochSeconds"]
        )
        t = s["timeSeconds"][keep]
        for field, color, style, label in (
            ("originalTorque", ".25", "-", "Recorded old torque"),
            ("A01OriginalPolicy", "#1469b0", "-", "A=0.1, original-equivalent rule"),
            (
                "A01FullFormulaCandidate",
                "#805098",
                "--",
                "A=0.1, full-formula candidate",
            ),
        ):
            axes[0, j].plot(
                t, s[field][keep], color=color, ls=style, lw=1.3, label=label
            )
        axes[0, j].set_title(f"{key[0]} | {key[1]} kg / {key[2]}")
        axes[0, j].set_ylabel("Net L5/S1 torque (N m)")
        for field, color, label in (
            ("oldPassiveBackTorque", "#287d4d", "Back passive"),
            ("oldPassiveAbdomenTorque", "#b66526", "Abdomen passive"),
            ("oldPassiveOtherTorque", ".4", "Other passive"),
        ):
            axes[1, j].plot(t, s[field][keep], color=color, lw=1.2, label=label)
        axes[1, j].set_ylabel("Recorded passive torque (N m)")
        for field, color, label in (
            ("A01FullBackTorque", "#287d4d", "Back"),
            ("A01FullAbdomenTorque", "#b66526", "Abdomen"),
            ("A01FullOtherTorque", ".4", "Other"),
        ):
            axes[2, j].plot(t, s[field][keep], color=color, lw=1.2, label=label)
        axes[2, j].set_ylabel("Full-candidate torque by group (N m)")
        axes[2, j].set_xlabel("Time from filtered-EMG start (s)")
        for i in range(3):
            axes[i, j].grid(alpha=0.15)
            axes[i, j].axhline(0, color=".6", ls=":", lw=0.8)
            axes[i, j].legend(fontsize=7, loc="best")
        data.append(
            dict(
                subject=row["subject"],
                trial=row["trial"],
                oldPeakNm=float(s["originalTorque"][keep].max()),
                A01OriginalPolicyPeakNm=float(s["A01OriginalPolicy"][keep].max()),
                A01FullFormulaPeakNm=float(s["A01FullFormulaCandidate"][keep].max()),
                A01FullFormulaMinimumNm=float(s["A01FullFormulaCandidate"][keep].min()),
                passiveMinimumNm=float(
                    sum(
                        s[k][keep]
                        for k in (
                            "oldPassiveBackTorque",
                            "oldPassiveAbdomenTorque",
                            "oldPassiveOtherTorque",
                        )
                    ).min()
                ),
            )
        )
    fig.suptitle(
        "CEINMS reference audit: changing A also changes torque; external calibration unavailable",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    fig.text(
        0.5,
        0.012,
        "Signed torques: extension positive. Native geometry and parameters retained. Purple is a model-structure diagnostic, not a corrected measured reference.",
        ha="center",
        fontsize=8,
    )
    fig.savefig(OUTPUT / "reference_torque_diagnostics.png", dpi=155)
    plt.close(fig)
    return data


def main():
    OUTPUT.mkdir(exist_ok=True, parents=True)
    DERIVED.mkdir(exist_ok=True, parents=True)
    review = json.loads(REVIEW.read_text())
    profile = json.loads(
        (ROOT / "data/processed/exp_activation_A01/profile.json").read_text()
    )
    configs = {
        r["subject"]: Path(r["originalConfiguration"]) for r in profile["profiles"]
    }
    records, signals = [], {}
    for row in review["basicRecords"] + review["supplementalNumberedRecords"]:
        if "error" in row:
            continue
        try:
            record, trace = replay(row, configs[row["subject"]])
        except Exception as error:
            record, trace = (
                dict(
                    subject=row["subject"],
                    trial=row["trial"],
                    candidateExported=False,
                    error=type(error).__name__ + ": " + str(error),
                ),
                None,
            )
        records.append(record)
        if trace is not None:
            signals[(row["subject"], row["trial"])] = trace
        print(
            row["subject"],
            row["trial"],
            "replayed" if trace is not None else record.get("error", "policy mismatch"),
            flush=True,
        )
    example = examples(review, records, signals)
    summary = []
    for subject in configs:
        rows = [r for r in records if r["subject"] == subject]
        good = [r for r in rows if r.get("candidateExported")]
        summary.append(
            dict(
                subject=subject,
                records=len(rows),
                verifiedOriginalEquivalentPolicy=len(good),
                maximumActiveReplayRMSE_N=max(
                    r["checks"]["reconstructedStaticActiveRMSE_N"] for r in good
                ),
                maximumElasticReplayRMSE_N=max(
                    r["checks"]["elasticMatchesPassiveExportRMSE_N"] for r in good
                ),
                minimumFullCandidateTorqueNm=min(
                    r["A01FullFormulaMinimumNm"] for r in good
                ),
                maximumFullCandidateTorqueNm=max(
                    r["A01FullFormulaPeakNm"] for r in good
                ),
            )
        )
    payload = dict(
        independentTorqueAvailable=False,
        independentlyCalibrated=False,
        originalDataModified=False,
        originalReferenceReplaced=False,
        userConfirmedNoIndependentForceOrIDData=True,
        candidateStatus="Sensitivity estimates on native exported geometry; not independently calibrated references.",
        stateRecovery="Back activation from recorded fibre geometry; abdominal activation from recorded active-force equation; this is NOT an independent forward validation.",
        unitPolicy="Replay the degree interpretation confirmed by recorded pennation; anatomical XML units need review before parameter recalibration.",
        fullFormulaVelocityPolicy="Recorded physical fibre velocity plus 15-point causal-spline derivative of the change in fibre length.",
        upstreamCoreCommit="5e53f3293ff4a5662ce575e83d2610c146de35ff",
        upstreamCoreURL="https://github.com/CEINMS-RT/ceinmsrt-core-cpp/blob/5e53f3293ff4a5662ce575e83d2610c146de35ff/lib/NMSmodel/MTU.cpp",
        summary=summary,
        examples=example,
        records=records,
    )
    (OUTPUT / "audit.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n"
    )
    fields = (
        "subject",
        "trial",
        "candidateExported",
        "oldPeakNm",
        "A01OriginalPolicyPeakNm",
        "A01FullFormulaPeakNm",
        "A01FullFormulaMinimumNm",
        "passiveTorqueMinimumNm",
        "passiveTorqueMaximumNm",
        "newActivationOutsidePhysicalRangeFraction",
        "error",
    )
    with (OUTPUT / "records.csv").open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: row.get(k) for k in fields} for row in records)
    print(json.dumps(dict(summary=summary, examples=example), indent=2), flush=True)


if __name__ == "__main__":
    main()
