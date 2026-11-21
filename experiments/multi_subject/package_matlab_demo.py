"""Package the frozen models and prepared recordings for the MATLAB demos.

This maintenance utility only copies data and extracts unchanged C equations.
The demos themselves need neither Python nor the experiment directories.
"""

from pathlib import Path
import hashlib
import json
import shutil
import sys

import numpy as np
from scipy.io import loadmat, savemat

from control_nine_cases import ROOT, BASE, load_subject, many, digest


def cells(values):
    array = np.empty((1, len(values)), dtype=object)
    for index, value in enumerate(values):
        array[0, index] = value
    return array


def matlab_values(value):
    """Keep MATLAB controller arithmetic in double precision, not int64."""
    if isinstance(value, dict):
        return {key: matlab_values(item) for key, item in value.items()}
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, list):
        if value and isinstance(value[0], dict):
            return cells([matlab_values(item) for item in value])
        return np.asarray(value, dtype=float)
    return value


def sequence_data(group):
    seq = group["sequence"]
    segments = []
    for segment in many(seq["segments"]):
        segments.append(
            {
                key: float(segment[key])
                for key in ("loadKg", "speedCode", "startIndex", "endIndex")
            }
        )
    return dict(
        time=np.asarray(seq["t"]).reshape(-1, 1),
        referenceTheta=np.asarray(seq["q"]).reshape(-1, 1),
        referenceTorque=np.asarray(seq["torque"]).reshape(-1, 1),
        torqueValid=np.asarray(
            seq.get("torqueReferenceValid", np.isfinite(seq["torque"])), dtype=bool
        ).reshape(-1, 1),
        referenceFilled=np.asarray(seq["referenceIsFilled"], dtype=bool).reshape(-1, 1),
        backActivation=np.asarray(seq["a"]),
        abdomenActivation=np.asarray(seq["front"]).reshape(-1, 1),
        mass=np.asarray(seq["m"]).reshape(-1, 1),
        nominalLoad=np.asarray(seq["loadKg"], dtype=float).reshape(-1, 1),
        waiting=np.asarray(seq["waiting"], dtype=float).reshape(-1, 1),
        dt=float(seq["dt"]),
        standingAngle=float(seq["standingAngle"]),
        initialState=np.array([seq["q"][0], seq["v0"]]),
        segments=cells(segments),
    )


def main():
    demo = ROOT / "Demo"
    data = demo / "BSE/data"
    control_data = demo / "NMPC/data"
    native = demo / "BSE/private/native"
    for folder in (data, native, control_data):
        folder.mkdir(parents=True, exist_ok=True)
    selected = json.loads((BASE / "selected_nmpc.json").read_text())
    manifest = dict(
        version=1,
        sourceReleaseSHA256=digest(BASE / "release_manifest.json"),
        subjects=[],
        independentValidation=False,
    )
    for subject, candidate in selected.items():
        package, models, groups, model_path, replay_path = load_subject(subject)
        all_groups = many(loadmat(replay_path, simplify_cells=True)["result"]["groups"])
        example = groups[0] if subject == "240528Pilot" else all_groups[-1]
        assert len(many(example["sequence"]["segments"])) == 3
        report_path = BASE / subject / "nmpc" / candidate / "evaluation.json"
        options = json.loads(report_path.read_text())["options"]
        clean_models = [
            dict(
                parameters=m["parameters"],
                neuralParameters=m["neuralParameters"],
                generator=m["generator"],
                intentDecoder="subtractive",
            )
            for m in models
        ]
        packed = dict(
            id=subject,
            shortName=subject[2:6][1:],
            balanceController=package["controller"],
            models=cells(clean_models),
            sequences=cells([sequence_data(g) for g in groups]),
            exampleSequence=sequence_data(example),
            modelSHA256=digest(model_path),
            sourceReplaySHA256=digest(replay_path),
            activationSources=cells(["Experimental activation", "MAP activation"]),
            mapProfilePolicy="Frozen outputs of the trained subject MAP; no EMG substitution or refit",
        )
        short = {"240528Pilot": "528", "240725Exp_A01": "725", "240730Exp_A01": "730"}[
            subject
        ]
        packed["shortName"] = short
        target = data / f"subject_{short}.mat"
        savemat(
            target, dict(subject=packed), do_compression=True, long_field_names=True
        )
        control_target = control_data / f"controller_{short}.mat"
        savemat(
            control_target,
            dict(controller=matlab_values(options), modelSHA256=digest(model_path)),
            do_compression=True,
            long_field_names=True,
        )
        manifest["subjects"].append(
            dict(
                subject=subject,
                file=target.relative_to(demo).as_posix(),
                fileSHA256=digest(target),
                modelSHA256=digest(model_path),
                controller=candidate,
                controllerFile=control_target.relative_to(demo).as_posix(),
                controllerFileSHA256=digest(control_target),
                conditions=9,
                fixedModelVectors=2,
                example=[
                    dict(loadKg=s["loadKg"], speedCode=s["speedCode"])
                    for s in many(example["sequence"]["segments"])
                ],
            )
        )

    shutil.copyfile(
        ROOT / "docs/nine_case_75pct/bse_parameter_values.csv",
        data / "bse_parameters.csv",
    )
    shutil.copyfile(
        ROOT / "docs/nine_case_75pct/subject_controller_parameters.csv",
        demo / "NMPC/controller_settings.csv",
    )

    if "--data-only" in sys.argv:
        manifest["nativeEquationOrigin"] = (
            "Unchanged frozen Hill, tendon, neural-control and joint equations; C wrappers reorganized only"
        )
        (demo / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        print("Updated the three subject data packages; native sources unchanged.")
        return

    # Keep only the shared equations; omit unused historical rollout gateways.
    old = ROOT / "experiments/hill_compliant"
    mechanics = (old / "hill_pair_core.c").read_text()
    mechanics = mechanics[: mechanics.index("/* The two muscle units")]
    (native / "hill_mechanics.h").write_text(mechanics)
    neural = (old / "hill_sequence_core.c").read_text()
    neural = neural[neural.index("/* Static force inversion") :]
    neural = neural[: neural.index("#if !defined(MATLAB_MEX_FILE)")]
    movement = (old / "hill_neural_core.c").read_text()
    movement = movement[movement.index("static void movement_command") :]
    movement = movement[: movement.index("#if !defined(MATLAB_MEX_FILE)")]
    decoder = "static double decode_intent_signal(double a, double lead, double tonic, double gain) {\n    return clamp01(a + lead - gain * tonic);\n}\n"
    (native / "hill_neural_control.h").write_text(neural + decoder + movement)
    plant = (ROOT / "experiments/legacy_mpc/hill_mpc_state.c").read_text()
    plant = plant[plant.index("#include <string.h>") :]
    plant = '#include "hill_mechanics.h"\n#include "hill_neural_control.h"\n' + plant
    (native / "hill_plant.h").write_text(plant)
    manifest["nativeEquationOrigin"] = (
        "Unchanged frozen Hill, tendon, neural-control and joint equations; C wrappers reorganized only"
    )
    (demo / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
