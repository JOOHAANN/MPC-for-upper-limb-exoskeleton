"""Install one current model/data package and controller per demo subject."""

import csv
import ctypes
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

import control_nine_cases as control
import recalibrate_bse_pearson as calibration
from activation_forecast import fit_predictor
from audit_subjects import ROOT
from control_nine_cases import digest, many, attach_planned_inputs
from native_library import load_library
from package_matlab_demo import sequence_data, matlab_values
from package_pearson_demo import reference_hash, restore_sequence
from prepare_subjects import structs
from run_hill_mpc import Plant as LegacyPlant

DEMO = ROOT / "Demo"
OUT = ROOT / "results/four_subject_demo"
SUBJECTS = {
    "528": "240528Pilot",
    "617": "240617Pilot",
    "725": "240725Exp_A01",
    "730": "240730Exp_A01",
}


class CurrentPlant(LegacyPlant):
    """Stateful access to the shared 32+16-parameter MATLAB/native plant."""

    def __init__(self, model, controller, angle):
        super().__init__(model, controller, angle)
        source = Path(__file__).with_name("demo_hill_state.c")
        dependencies = [source, *sorted((DEMO / "BSE/private/native").glob("*.h"))]
        library = load_library(source, dependencies, ROOT)
        for name in (
            "hill_mpc_init",
            "hill_mpc_advance",
            "hill_mpc_predict",
            "hill_mpc_predict_input_trend",
            "hill_mpc_predict_map",
        ):
            getattr(library, name).argtypes = getattr(self.library, name).argtypes
        self.library = library
        self.r = np.ascontiguousarray(np.r_[self.r, np.zeros(16 - len(self.r))])


def native_sequence(sequence):
    """Convert portable field names back to the calibration representation."""
    return dict(
        t=np.asarray(sequence["time"]),
        q=np.asarray(sequence["referenceTheta"]),
        torque=np.asarray(sequence["referenceTorque"]),
        a=sequence["backActivation"],
        front=sequence["abdomenActivation"],
        m=sequence["mass"],
        loadKg=sequence["nominalLoad"],
        waiting=sequence["waiting"],
        dt=sequence["dt"],
        v0=sequence["initialState"][1],
        standingAngle=sequence["standingAngle"],
        referenceIsFilled=sequence["referenceFilled"],
        torqueReferenceValid=sequence["torqueValid"],
        segments=sequence["segments"],
        stageId=np.zeros(len(sequence["time"])),
        speedCode=np.zeros(len(sequence["time"])),
    )


def main():
    control.Plant = CurrentPlant
    parameter_names = list(
        next(csv.DictReader((DEMO / "BSE/data/bse_parameters.csv").open()))
    )[3:]
    if len(parameter_names) == 42:
        parameter_names += [
            "frontIntentGain",
            "frontLoadIntentGain",
            "backLoadIntentGain",
            "backSquaredIntentGain",
            "loadDerivativeIntentGain",
            "intentBiasRad",
        ]
    parameter_rows = []
    fit_rows = []
    manifest = []
    for short, subject_id in SUBJECTS.items():
        folder = OUT / subject_id
        folder.mkdir(parents=True, exist_ok=True)
        if short == "617":
            source_file = folder / "calibration_input.mat"
            original = loadmat(source_file, simplify_cells=True)
            raw_groups = many(original["groups"])
            sequences = []
            for group in raw_groups[:9]:
                sequence = sequence_data(group)
                sequence["segments"] = structs(
                    [
                        dict(
                            segment,
                            **{
                                k: v
                                for k, v in raw.items()
                                if k in ("recordId", "cadenceBPM")
                            },
                        )
                        for segment, raw in zip(
                            many(sequence["segments"]),
                            many(group["sequence"]["segments"]),
                        )
                    ]
                )
                sequences.append(sequence)
            example = sequence_data(raw_groups[-1])
            example["segments"] = structs(
                [
                    dict(
                        segment,
                        **{
                            k: v
                            for k, v in raw.items()
                            if k in ("recordId", "cadenceBPM")
                        },
                    )
                    for segment, raw in zip(
                        many(example["segments"]),
                        many(raw_groups[-1]["sequence"]["segments"]),
                    )
                ]
            )
            subject = dict(
                id=subject_id,
                shortName=short,
                balanceController=original["controller"],
                sequences=structs(sequences),
                exampleSequence=example,
                mapModel=original["mapModel"],
                activationSources=structs(
                    ["Experimental activation", "MAP activation"]
                ),
                conditionPolicy="Nine complete recorded bouts; measured cadence; nominal speed unknown",
            )
            raw_groups = raw_groups[:9]
        else:
            original, _, raw_groups, _, replay_file = control.load_subject(subject_id)
            all_groups = many(
                loadmat(replay_file, simplify_cells=True)["result"]["groups"]
            )
            example_group = raw_groups[0] if short == "528" else all_groups[-1]
            sequences = [sequence_data(g) for g in raw_groups]
            example = sequence_data(example_group)
            subject = dict(
                id=subject_id,
                shortName=short,
                balanceController=original["controller"],
                activationSources=structs(
                    ["Experimental activation", "MAP activation"]
                ),
            )
        model_list = []
        passed = []
        for source in (1, 2):
            fitted = (
                ROOT / "results/pearson_recalibration" / subject_id / f"source{source}"
            )
            model = loadmat(fitted / "model.mat", simplify_cells=True)["model"]
            r = np.asarray(model["neuralParameters"]).reshape(-1)
            model = dict(
                parameters=model["parameters"],
                neuralParameters=np.r_[r, np.zeros(16 - len(r))],
                generator=model["generator"],
                intentDecoder="subtractive",
            )
            model_list.append(model)
            evaluation = json.loads((fitted / "evaluation.json").read_text())
            passed.append(evaluation["final"]["correlationTargetsPassed"])
            fit_rows.append(
                dict(Subject=subject_id, Source=source, **evaluation["final"])
            )
            parameter_rows.append(
                dict(
                    subject=subject_id,
                    source=source,
                    intentDecoder="subtractive",
                    **dict(
                        zip(
                            parameter_names,
                            np.r_[model["parameters"], model["neuralParameters"]],
                        )
                    ),
                )
            )
        model_file = folder / "current_models.mat"
        savemat(
            model_file,
            dict(models=structs(model_list), controller=subject["balanceController"]),
            do_compression=True,
            long_field_names=True,
        )
        subject["models"] = structs(model_list)
        subject["modelVersion"] = "current"
        subject["modelSHA256"] = digest(model_file)
        subject["calibrationTargetsPassed"] = np.asarray(passed, dtype=bool)
        subject["referenceSHA256"] = reference_hash(sequences)
        subject["sequences"] = structs([restore_sequence(s) for s in sequences])
        subject["exampleSequence"] = restore_sequence(example)
        target = DEMO / "BSE/data" / f"subject_{short}.mat"
        savemat(
            target, dict(subject=subject), do_compression=True, long_field_names=True
        )
        template_file = folder / "controller_template.json"
        if not template_file.exists():
            base_short = "528" if short == "617" else short
            template = loadmat(
                DEMO / "NMPC/data" / f"controller_{base_short}.mat", simplify_cells=True
            )["controller"]

            def plain(value):
                if isinstance(value, dict):
                    return {k: plain(v) for k, v in value.items()}
                if isinstance(value, np.ndarray):
                    return plain(value.tolist())
                if isinstance(value, (list, tuple)):
                    return [plain(v) for v in value]
                if isinstance(value, np.generic):
                    return value.item()
                return value

            template_file.write_text(json.dumps(plain(template), indent=2) + "\n")
        options = json.loads(template_file.read_text())
        for key in ("points", "blocks", "predictionSubsteps", "actualSubsteps"):
            options[key] = int(options[key])
        options["inputForecast"] = "arx"
        options["mapPhaseObserver"] = False
        options["angleGovernorUpperRad"] = 2.0
        signature = dict(
            modelSHA256=subject["modelSHA256"],
            referenceSHA256=subject["referenceSHA256"],
            actualSubsteps=10,
            points=options["points"],
            predictionDt=options["predictionDt"],
        )
        guard_file = folder / "prediction_guard.json"
        cached = json.loads(guard_file.read_text()) if guard_file.exists() else None
        for group in raw_groups:
            inputs = []
            predictions = []
            for source, model in enumerate(model_list):
                record = calibration.fit.record(group["sequence"], model["generator"])
                _, outputs = calibration.evaluate(
                    model["parameters"],
                    model["neuralParameters"],
                    subject["balanceController"],
                    [record],
                    source,
                    10,
                )
                out, bad, activation = outputs[0]
                assert not bad
                inputs.append(activation)
                predictions.append(dict(backTorque=out[:, 3]))
            group["inputs"] = np.asarray(inputs)
            group["predictions"] = structs(predictions)
        attach_planned_inputs(raw_groups, model_list)
        options["activationResidualPredictor"] = fit_predictor(
            raw_groups,
            options,
            signature,
            folder / "activation_residual_predictor.json",
        )
        if cached is None or cached["signature"] != signature:
            guard, relative, by_source = control.calibrate_guard(
                model_list, subject["balanceController"], raw_groups, options
            )
            cached = dict(
                signature=signature,
                guardNm=guard.tolist(),
                relativeGuard=relative.tolist(),
                bySource=by_source,
            )
            guard_file.write_text(json.dumps(cached, indent=2) + "\n")
        options["sourcePredictionModels"] = cached["bySource"]
        options["predictionErrorGuardNm"] = cached["guardNm"]
        options["predictionErrorGuardRelative"] = cached["relativeGuard"]
        equilibrium = []
        for model in model_list:
            _, first = CurrentPlant(model, subject["balanceController"], 1.5).init(
                np.zeros(2), 0.0, 1.0, [1.5, 0.0]
            )
            equilibrium.append(first[3])
        options["recoveryBackTargetNm"] = max(12.0, max(equilibrium) + 2.0)
        selected_file = folder / "selected_controller.json"
        if selected_file.exists():
            options.update(json.loads(selected_file.read_text()))
        savemat(
            DEMO / "NMPC/data" / f"controller_{short}.mat",
            dict(controller=matlab_values(options), modelSHA256=subject["modelSHA256"]),
            do_compression=True,
            long_field_names=True,
        )
        (folder / "current_controller.json").write_text(
            json.dumps(options, indent=2) + "\n"
        )
        manifest.append(
            dict(
                subject=subject_id,
                file=f"BSE/data/subject_{short}.mat",
                fileSHA256=digest(target),
                modelSHA256=subject["modelSHA256"],
                conditions=9,
                fixedParameterVectors=2,
                controllerFile=f"NMPC/data/controller_{short}.mat",
                controllerFileSHA256=digest(
                    DEMO / "NMPC/data" / f"controller_{short}.mat"
                ),
                correlationTargetsPassed=passed,
            )
        )
        print(
            subject_id,
            "packaged current BSE, trained activation forecaster and recalibrated prediction guards",
            flush=True,
        )
    for filename, rows in [
        ("bse_parameters.csv", parameter_rows),
        ("fit_summary.csv", fit_rows),
    ]:
        with (DEMO / "BSE/data" / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    (DEMO / "manifest.json").write_text(
        json.dumps(
            dict(
                version=3,
                modelVersion="current",
                subjects=manifest,
                originalDataChanged=False,
                independentValidation=False,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
