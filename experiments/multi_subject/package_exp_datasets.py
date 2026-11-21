"""Publish two derived A=0.1 datasets beside their original Data directories.

Only owned derived directories can be refreshed. Original subject files are
read-only; all copied and regenerated files receive provenance in the manifest.
Default Torque.sto retains the observed active-only policy and is explicitly a
model estimate. The complete-formula alternative is named as a candidate.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil

from audit_subjects import DATA, ROOT
from audit_ceinms_reference import replay
from review_exp_one_cycle import sha

SOURCE_INPUTS = {"emg.sto", "emgFilt.sto", "ik.sto", "lmt.sto", "ma_L5_S1_Flex_Ext.sto"}


def copy_file(source, target, records, role):
    target.parent.mkdir(parents=True, exist_ok=True)
    expected = sha(source)
    shutil.copy2(source, target)
    if sha(target) != expected or sha(source) != expected:
        raise ValueError("Copy integrity check failed: " + str(source))
    records.append(
        dict(source=str(source), destination=str(target), sha256=expected, role=role)
    )


def package(subject, activation_profile, review):
    source = (DATA / subject).resolve()
    target = (DATA / (subject + "_A01")).resolve()
    if target.parent != source.parent or target == source:
        raise ValueError("Derived folder must be a distinct sibling of the source")
    marker = target / "derived_dataset.json"
    if target.exists():
        if (
            not marker.is_file()
            or json.loads(marker.read_text())["sourceSubject"] != subject
        ):
            raise FileExistsError(
                "Refusing to overwrite an unowned folder: " + str(target)
            )
    else:
        target.mkdir()
    marker.write_text(
        json.dumps(
            dict(sourceSubject=subject, status="building", isDerivedDataset=True),
            indent=2,
        )
        + "\n"
    )
    profile = next(p for p in activation_profile["profiles"] if p["subject"] == subject)
    files, trials = [], []
    copy_file(
        Path(profile["configuration"]),
        target / "subject_A01.xml",
        files,
        "active_A01_subject_configuration",
    )
    for original in sorted(source.glob("*.xml")):
        copy_file(
            original,
            target / "Original_configuration" / original.name,
            files,
            "original_configuration_archive",
        )
    for original in sorted(source.glob("*.osim")):
        copy_file(
            original, target / original.name, files, "unchanged_source_geometry_model"
        )
    geometry = source / "Geometry"
    if geometry.exists():
        for original in sorted(geometry.rglob("*")):
            if original.is_file():
                copy_file(
                    original,
                    target / "Geometry" / original.relative_to(geometry),
                    files,
                    "unchanged_source_geometry_asset",
                )
    for row in review["basicRecords"] + review["supplementalNumberedRecords"]:
        if row["subject"] != subject or "error" in row:
            continue
        original_directory = source / row["trial"]
        new_directory = target / row["trial"]
        if not new_directory.resolve().is_relative_to(target):
            raise ValueError("Trial escapes the requested dataset directory")
        for original in sorted(original_directory.iterdir()):
            if not original.is_file() or original.suffix.lower() not in (
                ".sto",
                ".csv",
            ):
                continue
            destination = (
                new_directory / original.name
                if original.name in SOURCE_INPUTS
                else new_directory / "Original_CEINMS" / original.name
            )
            copy_file(
                original,
                destination,
                files,
                (
                    "unchanged_input"
                    if original.name in SOURCE_INPUTS
                    else "original_CEINMS_output_archive"
                ),
            )
        native = ROOT / "data/processed/exp_activation_A01" / subject / row["trial"]
        for name in ("Activations.sto", "sensor_activation_A01.sto"):
            copy_file(
                native / name, new_directory / name, files, "recomputed_A01_activation"
            )
        record, _ = replay(row, Path(profile["originalConfiguration"]), new_directory)
        if not record.get("candidateExported"):
            raise ValueError("Original-equivalent policy not verified: " + row["trial"])
        derived = Path(record["derivedDirectory"])
        for name in (
            "Torque_A01_original_policy.sto",
            "Torque_A01_plus_elastic_candidate.sto",
            "Torque_A01_full_formula_candidate.sto",
            "torque_sensitivity.mat",
        ):
            copy_file(
                derived / name,
                new_directory / name,
                files,
                "uncalibrated_torque_sensitivity_candidate",
            )
        trials.append(
            dict(
                trial=row["trial"],
                loadKg=row["loadKg"],
                speedCode=row["speedCode"],
                skeletalSamples=record["samples"],
                internalChecks=record["checks"],
                defaultTorquePolicy="A01_original_equivalent_active_only_model_estimate",
                independentlyCalibrated=False,
            )
        )
        print(subject, row["trial"], "packaged", flush=True)
    windows = ROOT / "data/processed/exp_activation_A01" / subject / "selected_windows"
    for original in sorted(windows.glob("*.mat")):
        copy_file(
            original,
            target / "selected_windows" / original.name,
            files,
            "reviewed_activation_window",
        )
    for original in (
        ROOT / "docs/exp_activation_A01" / (subject + "_A01_grid.png"),
        ROOT / "docs/exp_activation_A01/activation_before_after.png",
        ROOT / "docs/ceinms_reference_audit/reference_torque_diagnostics.png",
    ):
        copy_file(original, target / "reports" / original.name, files, "review_figure")
    readme = f"""# {subject}_A01：A=0.1 派生数据

本目录与原始 `{subject}` 并列，原始目录未修改。本次包含 {len(trials)} 份已处理记录；它们不是新增受试者。

`subject_A01.xml` 是本目录激活参数配置。原 XML 归档于 `Original_configuration/`。

每份 `CalibrationData/` 记录中：

- `Activations.sto`：从模型输入 EMG 重算的 A=0.1 激活。
- `sensor_activation_A01.sto`：完整 filtered EMG 的八通道估计，不能混作已验证的肌骨输出。
- `emg.sto`、`emgFilt.sto`、`ik.sto`、`lmt.sto` 和力臂：与原数据一致。
- `MuscleForces.sto`、主动/被动力、肌纤维及肌腱状态：在原生骨骼时刻按 A=0.1 回算，保留原导出所体现的简化规则。
- **`Torque.sto` 是 A=0.1、原等价主动项规则下的模型估计，未经独立实验校准。** 不把它称为实测真值。
- `skeletal_activation_A01.sto` 和 `skeletal_excitation_recovered.sto` 明确记录了异步肌骨时刻的状态恢复。
- `Torque_A01_plus_elastic_candidate.sto`、`Torque_A01_full_formula_candidate.sto` 分别用于被动项和完整公式敏感性检查，不能当作已经修正正确的 reference。
- 原 activation、Torque、肌肉力和状态位于同试次的 `Original_CEINMS/`，避免与新结果混用。

`selected_windows/` 保留前次审查的单动作窗口，`reports/` 为对照图。已有断档、负激活和部分肌骨导出缺失仍须复核；本目录没有补造缺失试验，也不等于 MAP/BSE/MPC 已重新训练。

原始输入的时间戳保留，新生成的肌骨状态使用原生骨骼时间戳；没有按试次重调 reference 角度。`derived_dataset.json` 列出逐文件来源、SHA-256 和求力规则。通用原始 subject 清单会跳过带此标记的派生目录，防止重复统计。

完整方法说明位于：

{ROOT / 'docs/CEINMS_REFERENCE_AUDIT.zh-CN.md'}

{ROOT / 'docs/EXP_ACTIVATION_A01.zh-CN.md'}
"""
    (target / "README.zh-CN.md").write_text(readme, encoding="utf-8")
    manifest = dict(
        sourceSubject=subject,
        sourceDirectory=str(source),
        directory=str(target),
        isDerivedDataset=True,
        status="complete",
        shapeFactor=0.1,
        includedRecords=len(trials),
        producer="experiments/multi_subject/package_exp_datasets.py",
        originalDataModified=False,
        independentlyCalibrated=False,
        referenceType="CEINMS model estimate under the original-equivalent active-only rule",
        recoveredStateWarning="Same-state sensitivity, not a replay of the unavailable historical DLL or independent force validation",
        copiedSources=files,
        trials=trials,
        fileSHA256={
            str(p.relative_to(target)): sha(p)
            for p in sorted(target.rglob("*"))
            if p.is_file() and p != marker
        },
    )
    marker.write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print("COMPLETE", target, len(trials), flush=True)
    return dict(
        sourceSubject=subject,
        directory=str(target),
        records=len(trials),
        files=len(manifest["fileSHA256"]),
        manifestSHA256=sha(marker),
    )


def main():
    profile = json.loads(
        (ROOT / "data/processed/exp_activation_A01/profile.json").read_text()
    )
    review = json.loads((ROOT / "docs/exp_one_cycle_review/review.json").read_text())
    summary = [
        package(subject, profile, review) for subject in ("240725Exp", "240730Exp")
    ]
    (ROOT / "docs/ceinms_reference_audit/dataset_locations.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
