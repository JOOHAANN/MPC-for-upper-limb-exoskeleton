"""Publish only 5/10/15 kg BSE/NMPC scenarios, retaining five-load MAP figures."""

from __future__ import annotations

import json
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT
from fit_subject import many
from quality_gate import sha
from run_exp_subjects import status
from report_subjects import sequence_figure, case_figure, summary_rows

PUBLISH = ROOT / "docs/exp_three_loads"
LOADS = (5, 10, 15)


def calibration_figure(d):
    fig, axes = plt.subplots(4, 3, figsize=(13.5, 10), sharex="col")
    for col, load in enumerate(LOADS):
        g, seg = next(
            (g, s)
            for g in many(d["r"]["groups"])[:-1]
            for s in many(g["sequence"]["segments"])
            if s["loadKg"] == load and s["split"] == "train"
        )
        seq = g["sequence"]
        ids = slice(int(seg["startIndex"]) - 1, int(seg["endIndex"]))
        t = seq["t"][ids] - seq["t"][ids.start]
        for source in (0, 1):
            p = many(g["predictions"])[source]
            nominal = many(g["nominal"])[source]
            for offset, key in ((0, "q"), (1, "torque")):
                ax = axes[2 * source + offset, col]
                ax.plot(
                    t,
                    seq[key][ids],
                    color="#cc3434",
                    ls="--",
                    lw=1.2,
                    label="Recorded IK / available CEINMS estimate",
                )
                ax.plot(
                    t,
                    nominal[key][ids],
                    color="#c8841c",
                    lw=1.5,
                    label="Fitted BSE, no assistance",
                )
                ax.plot(
                    t, p[key][ids], color="#1469b0", lw=1.2, label="Same BSE + NMPC"
                )
                ax.set_xlim(0, t[-1])
                ax.grid(alpha=0.2)
                if key == "q":
                    ax.set_ylim(-0.2, 1.85)
                if col == 0:
                    ax.set_ylabel(
                        ("Experimental" if source == 0 else "MAP")
                        + (" theta (rad)" if key == "q" else " net torque (N m)")
                    )
        axes[0, col].set_title(
            f"{load} kg / speed code {seg['speedCode']}", fontsize=11
        )
        axes[-1, col].set_xlabel("Action time (s)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8)
    fig.suptitle(
        d["subject"]
        + " | requested 5 / 10 / 15 kg calibration actions\nFitted actions; additional BSE holdouts are reported separately",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    return fig


def item(subject, folder):
    replay = loadmat(folder / "hill_mpc.mat", simplify_cells=True)["result"]
    groups = [
        g
        for g in many(replay["groups"])
        if all(s["loadKg"] in LOADS for s in many(g["sequence"]["segments"]))
    ]
    replay["groups"] = groups
    prep = json.loads((folder / "preparation.json").read_text())
    prep["cycles"] = [r for r in prep["cycles"] if r["loadKg"] in LOADS]
    meta = json.loads((folder / "mpc_evaluation.json").read_text())
    meta["metrics"] = [r for r in meta["metrics"] if r["loadKg"] in LOADS]
    return dict(
        subject=subject,
        r=replay,
        cache=loadmat(folder / "neural_models.mat", simplify_cells=True),
        meta=meta,
        prep=prep,
        bse=[
            r
            for r in json.loads((folder / "bse_metrics.json").read_text())
            if r["loadKg"] in LOADS
        ],
        policy=json.loads((folder / "fit_policy.json").read_text()),
        validation=json.loads((folder / "validation.json").read_text()),
        showUnassistedBoundary=True,
    )


def main():
    PUBLISH.mkdir(parents=True, exist_ok=True)
    folders = {
        "240725Exp_A01": ROOT / "results/exp_three_loads_nmpc/240725Exp_A01",
        "240730Exp_A01": ROOT / "results/exp_subjects_A01_nmpc/240730Exp_A01",
    }
    data = [item(subject, folder) for subject, folder in folders.items()]
    rows = summary_rows(data)
    portable = ROOT / "data/processed/exp_subjects_A01"
    for subject, folder in folders.items():
        shutil.copy2(folder / "neural_models.mat", portable / (subject + "_models.mat"))
        for name in (
            "fit_policy.json",
            "preparation.json",
            "bse_metrics.json",
            "mpc_evaluation.json",
            "validation.json",
        ):
            shutil.copy2(folder / name, PUBLISH / (subject + "_" + name))
        for suffix in (".png", ".json"):
            shutil.copy2(
                ROOT / "results/exp_subjects_A01/map" / (subject + "_map" + suffix),
                PUBLISH / (subject + "_map" + suffix),
            )
    shutil.copy2(
        ROOT / "results/exp_three_loads/three_anchor_inputs.png",
        PUBLISH / "three_anchor_inputs.png",
    )
    shutil.copy2(
        ROOT / "results/exp_three_loads/quality_gate.json",
        PUBLISH / "three_anchor_data_policy.json",
    )
    shutil.copy2(
        ROOT / "results/exp_three_loads/preserved_models.json",
        PUBLISH / "preserved_models.json",
    )
    with PdfPages(PUBLISH / "three_load_results.pdf") as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
        for source, color, offset in (
            ("Experimental activation", "#1469b0", -0.16),
            ("MAP activation", "#23844a", 0.16),
        ):
            chosen = [
                next(
                    r
                    for r in rows
                    if r["subject"] == d["subject"]
                    and r["source"] == source
                    and r["split"] == "train"
                )
                for d in data
            ]
            axes[0].bar(
                np.arange(2) + offset,
                [r["mpcThetaCorrelation"] for r in chosen],
                width=0.3,
                color=color,
                label=source,
            )
            axes[1].bar(
                np.arange(2) + offset,
                [r["worstBackExcessNm"] for r in chosen],
                width=0.3,
                color=color,
            )
        for ax in axes:
            ax.set_xticks(range(2), ["725 Exp", "730 Exp"])
            ax.grid(axis="y", alpha=0.2)
        axes[0].set(
            ylabel="Theta trend correlation",
            title="Three BSE calibration anchors: 5 / 10 / 15 kg",
        )
        axes[1].set(
            ylabel="Maximum excess back torque (N m)",
            title="NMPC calibration replay; extra actions reported separately",
        )
        axes[0].legend(fontsize=8)
        fig.suptitle(
            "MAP trained from excitation | BSE/NMPC requested loads only | completed subjects preserved"
        )
        fig.tight_layout()
        fig.savefig(PUBLISH / "completion_overview.png", dpi=160)
        pdf.savefig(fig)
        plt.close(fig)
        for d in data:
            if d["subject"] == "240725Exp_A01":
                calibration = calibration_figure(d)
                calibration.savefig(
                    PUBLISH / (d["subject"] + "_calibration.png"),
                    dpi=155,
                    bbox_inches="tight",
                )
                pdf.savefig(calibration, bbox_inches="tight")
                plt.close(calibration)
            for fig, slug in (
                (sequence_figure(d), "continuous"),
                (case_figure(d, 0), "cases_experimental"),
                (case_figure(d, 1), "cases_map"),
            ):
                fig.savefig(
                    PUBLISH / (d["subject"] + "_" + slug + ".png"),
                    dpi=155,
                    bbox_inches="tight",
                )
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
    result = dict(
        subjects=status(),
        BSELoadsKg=list(LOADS),
        MAPLoadsKg=[5, 7.5, 10, 12.5, 15],
        MAPModelsRefitted=False,
        threeLoadAnchorPolicy="One complete action for each requested load; additional actions are optional BSE holdouts.",
        missingTorqueSamplesUsedForFitting=False,
        referenceTorqueIndependentlyCalibrated=False,
        summary=rows,
        sources={s: str(p) for s, p in folders.items()},
        MAPHoldoutsAndBSEHoldoutsAreDifferent=True,
        notes=[
            "725 extra actions were excluded from BSE fitting but some were used by the already-frozen MAP; do not call these end-to-end pipeline holdouts.",
            "730 historical 12.5 kg replay is preserved on disk but excluded from this requested three-load report.",
        ],
    )
    (PUBLISH / "completion_status.json").write_text(json.dumps(result, indent=2) + "\n")
    (PUBLISH / "export_metadata.json").write_text(
        json.dumps(
            dict(
                files={
                    p.name: sha(p)
                    for p in sorted(PUBLISH.iterdir())
                    if p.is_file() and p.name != "export_metadata.json"
                }
            ),
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == "__main__":
    main()
