"""Package Pearson-calibrated models and explicitly identified unmet targets.

Original datasets, MAP predictions, references, released parameter packages
and existing NMPC settings remain intact. A completion marker is written last.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

from audit_subjects import ROOT
from control_nine_cases import digest, many
from prepare_subjects import structs
from recalibrate_bse_pearson import OUT

SUBJECTS = ("240528Pilot", "240725Exp_A01", "240730Exp_A01")
SOURCES = ("Experimental activation", "MAP activation")
DEMO = ROOT / "Demo"


def restore_sequence(sequence):
    sequence = dict(sequence)
    for name in (
        "time",
        "referenceTheta",
        "referenceTorque",
        "torqueValid",
        "referenceFilled",
        "abdomenActivation",
        "mass",
        "nominalLoad",
        "waiting",
    ):
        sequence[name] = np.asarray(sequence[name]).reshape(-1, 1)
    sequence["segments"] = structs(many(sequence["segments"]))
    return sequence


def reference_hash(sequences):
    checksum = hashlib.sha256()
    for sequence in sequences:
        for name in (
            "time",
            "referenceTheta",
            "referenceTorque",
            "torqueValid",
            "backActivation",
            "abdomenActivation",
            "mass",
            "waiting",
        ):
            checksum.update(name.encode())
            checksum.update(np.asarray(sequence[name], dtype="<f8").tobytes())
    return checksum.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include-unmet-730-experimental",
        action="store_true",
        help="Retain the nonpassing 730 experimental candidate for review; do not mark it passed.",
    )
    args = parser.parse_args()
    bundles = []
    summary_rows, condition_rows, parameter_rows = [], [], []
    with (DEMO / "BSE/data/bse_parameters.csv").open(encoding="utf-8-sig") as stream:
        parameter_names = list(next(csv.DictReader(stream)))[3:]
    assert len(parameter_names) == 42
    parameter_names += [
        "frontIntentGain",
        "frontLoadIntentGain",
        "backLoadIntentGain",
        "backSquaredIntentGain",
        "loadDerivativeIntentGain",
        "intentBiasRad",
    ]
    for subject_id in SUBJECTS:
        short = subject_id[3:6]
        original_file = DEMO / "BSE/data" / f"subject_{short}.mat"
        subject = loadmat(original_file, simplify_cells=True)["subject"]
        sequences = many(subject["sequences"])
        original_reference_hash = reference_hash(sequences)
        models = []
        source_passed = []
        for source in range(2):
            folder = OUT / subject_id / f"source{source+1}"
            evaluation = json.loads((folder / "evaluation.json").read_text())
            passed = bool(evaluation["final"]["passed"])
            allowed_unmet = (
                args.include_unmet_730_experimental
                and subject_id == "240730Exp_A01"
                and source == 0
            )
            assert (
                passed or allowed_unmet
            ), f"Unfinished fit: {subject_id} source {source+1}"
            source_passed.append(passed)
            assert len(evaluation["rows"]) == 9
            model = loadmat(folder / "model.mat", simplify_cells=True)["model"]
            model["intentDecoder"] = "subtractive"
            neural = np.asarray(model["neuralParameters"]).reshape(-1)
            model["neuralParameters"] = np.r_[neural, np.zeros(16 - len(neural))]
            models.append(model)
            before = evaluation["baseline"]
            after = evaluation["final"]
            summary_rows.append(
                dict(
                    Subject=subject_id,
                    Source=SOURCES[source],
                    TorqueMeanRBefore=before["meanTorqueR"],
                    TorqueMeanRAfter=after["meanTorqueR"],
                    ThetaMeanRBefore=before["meanThetaR"],
                    ThetaMeanRAfter=after["meanThetaR"],
                    TorqueRMSEBeforeNm=before["meanTorqueRMSE"],
                    TorqueRMSEAfterNm=after["meanTorqueRMSE"],
                    TorqueMeanRRequired=evaluation["targetTorqueMeanR"],
                    ThetaMeanRRequired=evaluation["targetThetaMeanR"],
                    MinimumThetaRad=after["minimumTheta"],
                    MaximumThetaRad=after["maximumTheta"],
                    MaximumCaseMeanAbdomenBackRatio=after["maximumAbdomenBackRatio"],
                    BoundaryPeakNm=after["boundaryPeakNm"],
                    ConditionCount=9,
                    Passed=passed,
                    Status=(
                        "Targets met"
                        if passed
                        else "Targets unmet; retained for review"
                    ),
                )
            )
            for row in evaluation["rows"]:
                condition_rows.append(
                    dict(Subject=subject_id, Source=SOURCES[source], **row)
                )
            vector = np.r_[model["parameters"], model["neuralParameters"]]
            parameter_rows.append(
                dict(
                    subject=subject_id,
                    source=SOURCES[source],
                    intentDecoder="subtractive",
                    **dict(zip(parameter_names, vector)),
                )
            )
        combined_file = OUT / subject_id / "neural_models.mat"
        savemat(
            combined_file,
            dict(
                subject=subject_id,
                models=structs(models),
                controller=subject["balanceController"],
                modelVersion="pearson",
                calibrationConditions=9,
                referenceSHA256=original_reference_hash,
            ),
            do_compression=True,
            long_field_names=True,
        )
        subject["models"] = structs(models)
        subject["originalModelSHA256"] = subject["modelSHA256"]
        subject["modelSHA256"] = digest(combined_file)
        subject["modelVersion"] = "pearson"
        subject["calibrationTargetsPassed"] = np.asarray(source_passed, dtype=bool)
        subject["referenceSHA256"] = original_reference_hash
        subject["sequences"] = structs([restore_sequence(s) for s in sequences])
        subject["exampleSequence"] = restore_sequence(subject["exampleSequence"])
        assert reference_hash(sequences) == original_reference_hash
        bundles.append(
            (
                short,
                subject,
                dict(
                    subject=subject_id,
                    originalPackageSHA256=digest(original_file),
                    modelSHA256=digest(combined_file),
                    referenceSHA256=original_reference_hash,
                    oneVectorPerSourceAcrossNineCases=True,
                    sourceTargetsPassed=source_passed,
                ),
            )
        )
    folder = DEMO / "BSE/data/pearson"
    folder.mkdir(parents=True, exist_ok=True)
    artifacts = []
    for short, subject, metadata in bundles:
        target = folder / f"subject_{short}.mat"
        savemat(
            target, dict(subject=subject), do_compression=True, long_field_names=True
        )
        artifacts.append(dict(**metadata, file=target.name, fileSHA256=digest(target)))
    for filename, rows in [
        ("fit_summary.csv", summary_rows),
        ("condition_metrics.csv", condition_rows),
        ("bse_parameters.csv", parameter_rows),
    ]:
        with (folder / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    ready = dict(
        ready=True,
        passed=all(row["Passed"] for row in summary_rows),
        passingSourceModels=sum(row["Passed"] for row in summary_rows),
        unmetTargetsRetainedForReview=[
            dict(subject=row["Subject"], source=row["Source"])
            for row in summary_rows
            if not row["Passed"]
        ],
        version="pearson",
        subjects=artifacts,
        criteria="Unweighted mean of nine case-wise Pearson r values; finite reference support only",
        referenceChanged=False,
        MAPProfilesChanged=False,
        HillMechanicsChanged=False,
        neuralPolicyExtended=True,
        independentValidation=False,
        NMPCStillUsesReleasedModels=True,
        optimizationCodeSHA256=digest(
            Path(__file__).with_name("recalibrate_bse_pearson.py")
        ),
        summary=summary_rows,
    )
    (folder / "READY.json").write_text(json.dumps(ready, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                ready=True,
                passed=ready["passed"],
                passingSourceModels=ready["passingSourceModels"],
                subjects=3,
                sourceModels=6,
                comparisons=54,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
