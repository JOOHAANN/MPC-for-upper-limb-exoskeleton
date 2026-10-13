"""Publish verified scientific figures and numerical tables into docs/.

Run after demo_hill_compliant and test_hill_compliant have completed. The
reference-invariance check prevents exporting a rebuilt/rescaled reference.
"""
import csv
import shutil
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from package_hill_results import ROOT

OUT = ROOT / "docs/hill_compliant_evaluation"
RESULTS = ROOT / "results/hill_compliant"


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(name, rows):
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def number(row, key):
    return float(row[key])


def summarize(mode, records):
    summaries = []
    for source in ("Experimental activation", "MAP activation"):
        rows = [row for row in records if row["source"] == source]
        assert len(rows) == 9
        values = lambda key: np.array([number(row, key) for row in rows])
        timing = values("valleyTimeMAE")
        summaries.append({
            "mode": mode, "source": source, "trials": len(rows),
            "thetaRMSEDegrees": values("thetaRMSE").mean()*180/np.pi,
            "torqueRMSENm": values("torqueRMSE").mean(),
            "meanTorqueCorrelation": values("torqueCorrelation").mean(),
            "meanCorrelation": values("correlation").mean(),
            "twoValleyMatches": int(values("Wmatched").sum()),
            "valleyTimeMAESecondsMatchedOnly": np.nanmean(timing),
            "minimumDegrees": values("minAngle").min()*180/np.pi,
            "maximumDegrees": values("maxAngle").max()*180/np.pi,
            "backPassivePeakNm": values("backPassiveMax").max(),
            "frontPassivePeakNm": values("frontPassiveMax").max(),
            "boundaryPeakNm": values("boundaryMax").max(),
            "boundaryRMSMeanNm": values("boundaryRMS").mean(),
            "gravityPeakNm": values("gravityMax").max(),
            "jointDampingPeakNm": values("dampingMax").max(),
            "outsideCheckRangeMeanPercent": values("outsidePercent").mean(),
            "forceResidualMaximum": values("equilibriumResidualMax").max(),
            "momentumResidualMaximumNm": values("momentumResidualMax").max(),
        })
    return summaries


def parameter_rows(package):
    # This table exposes the actual calibration rather than generic defaults.
    from fit_hill_compliant import PARAMETER_NAMES
    records = []
    parameters = lambda model: np.asarray(model["parameters"] if isinstance(model,dict) else model.parameters).reshape(-1)
    for source, model in enumerate(package["models"], 1):
        records.append({"mode": "Calibration", "source": source, "heldSpeed": 0,
                        **dict(zip(PARAMETER_NAMES, parameters(model)))})
    for source in (1, 2):
        for fold, speed in enumerate((35, 45, 55)):
            model = package["validationModels"][source-1, fold]
            records.append({"mode": "Whole-speed holdout", "source": source, "heldSpeed": speed,
                            **dict(zip(PARAMETER_NAMES, parameters(model)))})
    return records


def main():
    package = loadmat(ROOT / "data/processed/hill_compliant/calibrated_models.mat", simplify_cells=True)
    legacy = loadmat(ROOT / "data/processed/hill_forward/calibrated_models.mat", simplify_cells=True)
    for current, original in zip(package["trials"], legacy["trials"]):
        assert np.array_equal(current["q"], original["q"])
        assert current["v0"] == original["v0"]
    from prepare_hill_trials import GROUPS, WINDOWS
    originals = np.vstack([loadmat(ROOT / "data/processed" /
        f"training_map__muscle_tendon__{group}_done__build_data/Combination_raw.mat")["Combination"] for group in GROUPS])
    audit = []
    for current, (left, right) in zip(package["trials"], WINDOWS):
        selected = originals[np.arange(left, right-1, 10), 2]
        assert np.array_equal(current["q"], selected)
        audit.append({"loadKg": current["loadKg"], "speed": current["speed"],
            "samples": len(selected), "exactOriginalReference": True,
            "maximumReferenceDifferenceRad": float(np.max(np.abs(current["q"]-selected)))})
    convergence = read_csv(RESULTS / "convergence.csv")
    assert len(convergence) == 36
    assert max(number(row, "angleDifference") for row in convergence) < .01
    assert max(number(row, "torqueDifference") for row in convergence) < 3
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv("reference_audit.csv", audit)
    summary = []
    for mode, folder, name in (
        ("Calibration", RESULTS, "calibration_metrics.csv"),
        ("Whole-speed holdout", RESULTS / "validation", "holdout_metrics.csv"),
    ):
        metrics = read_csv(folder / "metrics.csv")
        summary.extend(summarize(mode, metrics))
        shutil.copy2(folder / "metrics.csv", OUT / name)
    write_csv("summary.csv", summary)
    write_csv("parameters.csv", parameter_rows(package))
    for name in ("convergence.csv", "stability.csv"):
        shutil.copy2(RESULTS / name, OUT / name)
    for speed in (35, 45, 55):
        for name in (f"comparison_speed{speed}", f"theta_focus_speed{speed}", f"soft_limits_speed{speed}"):
            for suffix in (".png", ".pdf"):
                shutil.copy2(RESULTS / "validation" / (name+suffix), OUT / (name+suffix))
        shutil.copy2(RESULTS / f"theta_focus_speed{speed}.png", OUT / f"calibration_theta_speed{speed}.png")
    for name in ("torque_components_speed45.png", "activation_calibration_speed45.png", "equation_budget_speed45.png"):
        shutil.copy2(RESULTS / "validation" / name, OUT / name)
    from fit_hill_compliant import LOW, HIGH, PARAMETER_NAMES, SHARED
    units = ["Nm", "1", "1", "Nm s/rad", "L0/s", "1", "1", "s", "Nm", "1", "1", "1",
             "m", "rad", "1", "1", "1", "1", "kg", "m", "m", "m", "m", "m", "L0/s", "1", "s", "rad",
             "Nm/rad", "Nm", "rad", "rad"]
    write_csv("parameter_bounds.csv", [{"parameter": name, "unit": units[i], "lower": LOW[i],
               "upper": HIGH[i], "sharedBetweenSources": i in SHARED} for i, name in enumerate(PARAMETER_NAMES)])
    details = []
    for record in parameter_rows(package):
        for i, name in enumerate(PARAMETER_NAMES):
            value = float(record[name])
            details.append({"mode": record["mode"], "source": record["source"], "heldSpeed": record["heldSpeed"],
                            "parameter": name, "value": value, "unit": units[i],
                            "atBound": min(value-LOW[i], HIGH[i]-value)/(HIGH[i]-LOW[i]) < .001})
    write_csv("parameter_details.csv", details)
    derived = []
    for record in parameter_rows(package):
        for unit in ("back", "front"):
            arm = float(record[f"{unit}MomentArm"])
            optimal = arm / float(record[f"{unit}ArmToLength"])
            ref = float(record[f"{unit}ReferenceLength"])
            slack = float(record[f"{unit}ReferenceMTLength"]) - optimal * np.sqrt(ref**2-np.sin(float(record["referencePennation"]))**2)
            derived.append({"mode": record["mode"], "source": record["source"], "heldSpeed": record["heldSpeed"],
                "unit": unit, "momentArmM": arm, "optimalFiberLengthM": optimal,
                "tendonSlackLengthM": slack, "maximumIsometricForceN": float(record[f"{unit}Capacity"])/arm})
    write_csv("derived_muscle_parameters.csv", derived)
    before = ROOT / "results/hill_geometry_baseline"
    comparison = []
    for mode, filename, current_path in (
        ("Calibration", "calibration_metrics.csv", RESULTS / "metrics.csv"),
        ("Whole-speed holdout", "holdout_metrics.csv", RESULTS / "validation/metrics.csv"),
    ):
        if not (before / filename).exists():
            continue
        for version, rows in (("Previous fixed M/L1", read_csv(before / filename)),
                              ("Optimized geometry and soft limits", read_csv(current_path))):
            for source in ("Experimental activation", "MAP activation"):
                selected = [row for row in rows if row["source"] == source]
                comparison.append({"mode": mode, "version": version, "source": source,
                    "thetaRMSEDegrees": np.mean([number(row,"thetaRMSE") for row in selected])*180/np.pi,
                    "torqueRMSENm": np.mean([number(row,"torqueRMSE") for row in selected]),
                    "twoValleyMatches": int(sum(number(row,"Wmatched") for row in selected)),
                    "boundaryPeakNm": max(number(row,"boundaryMax") for row in selected)})
    if comparison:
        write_csv("before_after.csv", comparison)
    for row in summary:
        print(row)
    print("Verified figures and tables:", OUT)


if __name__ == "__main__":
    main()
