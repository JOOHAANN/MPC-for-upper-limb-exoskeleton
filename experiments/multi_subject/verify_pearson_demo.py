"""Check the portable MATLAB replay against the calibration and frozen data."""

import csv
import json

import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT
from control_nine_cases import digest, many
from package_pearson_demo import reference_hash


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    demo = ROOT / "Demo"
    folder = demo / "BSE/data/pearson"
    ready = json.loads((folder / "READY.json").read_text())
    matlab = read_csv(demo / "BSE/output/pearson/bse_metrics.csv")
    fitted = read_csv(folder / "condition_metrics.csv")

    def key(row):
        return (
            row["Subject"],
            row["Source"],
            float(row.get("LoadKg", row.get("loadKg"))),
            float(row.get("SpeedCode", row.get("speedCode"))),
        )

    lookup = {key(row): row for row in fitted}
    assert len(matlab) == len(lookup) == 54
    differences = []
    for row in matlab:
        other = lookup[key(row)]
        differences.append(
            [
                abs(float(row["TorqueCorrelation"]) - float(other["torqueR"])),
                abs(float(row["ThetaCorrelation"]) - float(other["thetaR"])),
                abs(float(row["TorqueRMSENm"]) - float(other["torqueRMSE"])),
            ]
        )
    maximum = np.max(differences, axis=0)
    assert maximum[0] < 1e-6 and maximum[1] < 1e-6 and maximum[2] < 1e-4

    for short in ("528", "725", "730"):
        original = loadmat(demo / f"BSE/data/subject_{short}.mat", simplify_cells=True)[
            "subject"
        ]
        updated = loadmat(folder / f"subject_{short}.mat", simplify_cells=True)[
            "subject"
        ]
        assert reference_hash(many(original["sequences"])) == reference_hash(
            many(updated["sequences"])
        )
        assert reference_hash([original["exampleSequence"]]) == reference_hash(
            [updated["exampleSequence"]]
        )
        models = many(updated["models"])
        assert len(models) == 2
        for model in models:
            assert np.asarray(model["parameters"]).size == 32
            assert np.asarray(model["neuralParameters"]).size == 16

    release = json.loads(
        (ROOT / "results/nine_case_75pct/release_manifest.json").read_text()
    )
    frozen = dict(release["artifactSHA256"], **release["sourceCodeSHA256"])
    for name, expected in frozen.items():
        assert digest(ROOT / name) == expected, f"Frozen file changed: {name}"

    diagnostics = read_csv(demo / "BSE/diagnostics/pearson/muscle_boundary_metrics.csv")
    assert len(diagnostics) == 54
    report = dict(
        replayChecksPassed=True,
        allCorrelationTargetsPassed=ready["passed"],
        passingSourceModels=ready["passingSourceModels"],
        unmetTargetsRetainedForReview=ready["unmetTargetsRetainedForReview"],
        conditionSourceComparisons=54,
        maximumTorqueCorrelationDifference=float(maximum[0]),
        maximumThetaCorrelationDifference=float(maximum[1]),
        maximumTorqueRMSEDifferenceNm=float(maximum[2]),
        frozenFilesUnchanged=len(frozen),
        referencesAndPreparedActivationsUnchanged=True,
        sourceSpecificFixedParameterVectors=6,
        independentValidation=False,
        nmpcModelVersion="frozen",
        nmpcRevalidatedOnPearsonModels=False,
    )
    (folder / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    verification = json.loads((demo / "verification.json").read_text())
    verification["scope"] = "Previously released frozen BSE and NMPC numerical checks"
    verification["pearsonRecalibration"] = report
    (demo / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    manifest = json.loads((demo / "manifest.json").read_text())
    manifest["version"] = 2
    manifest["nativeEquationOrigin"] = (
        "Frozen Hill mechanics retained. Optional neural intent terms are enabled by "
        "Pearson calibration parameters and are zero for the frozen NMPC models."
    )
    manifest["pearsonRecalibration"] = dict(
        marker="BSE/data/pearson/READY.json",
        markerSHA256=digest(folder / "READY.json"),
        verification="BSE/data/pearson/verification.json",
        allCorrelationTargetsPassed=ready["passed"],
        passingSourceModels=ready["passingSourceModels"],
        nmpcModelVersion="frozen",
    )
    for name in manifest["runtimeSourceSHA256"]:
        manifest["runtimeSourceSHA256"][name] = digest(demo / name)
    manifest["runtimeSourceSHA256"]["BSE/inspect_bse_model.m"] = digest(
        demo / "BSE/inspect_bse_model.m"
    )
    (demo / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
