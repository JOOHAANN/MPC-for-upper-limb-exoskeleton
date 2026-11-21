"""Compare saved controller versions, including every post-action transient."""

import csv
import json
import shutil
import numpy as np
from scipy.io import loadmat
from matplotlib.backends.backend_pdf import PdfPages

from plot_hill_mpc_results import plt, ROOT, BLUE, AMBER, GREEN, RED


def main():
    before = loadmat(
        ROOT / "results/hill_mpc_before_recovery/hill_mpc.mat", simplify_cells=True
    )["result"]
    after = loadmat(ROOT / "results/hill_mpc/hill_mpc.mat", simplify_cells=True)[
        "result"
    ]
    out = ROOT / "results/hill_mpc"
    target = ROOT / "docs/hill_mpc_evaluation"
    rows = []
    figures = []
    for ga, gb in zip(before["groups"], after["groups"]):
        seq = gb["sequence"]
        n = len(seq["t"])
        assert np.array_equal(seq["q"], ga["sequence"]["q"])
        for source in (0, 1):
            old, new = ga["predictions"][source], gb["predictions"][source]
            for i, seg in enumerate(seq["segments"]):
                j = int(seg["endIndex"])
                stop = min(n, j + 100)
                ids = slice(j, stop)
                rows.append(
                    {
                        "group": int(gb["group"]),
                        "source": source + 1,
                        "loadKg": int(seg["loadKg"]),
                        "speedCode": int(seg["speedCode"]),
                        "movementEndSeconds": float(seq["t"][j - 1]),
                        "beforeBackPeakNm": float(old["backTorque"][ids].max()),
                        "afterBackPeakNm": float(new["backTorque"][ids].max()),
                        "beforeRecoveryAngleRMSE": float(
                            np.sqrt(np.mean((old["q"][ids] - 1.5) ** 2))
                        ),
                        "afterRecoveryAngleRMSE": float(
                            np.sqrt(np.mean((new["q"][ids] - 1.5) ** 2))
                        ),
                        "afterRecoveryAssistPeakNm": float(
                            new["assistanceTorque"][ids].max()
                        ),
                        "afterAbdomenPeakNm": float(new["frontTorque"][ids].max()),
                        "afterBoundaryPeakNm": float(
                            np.abs(new["boundaryTorque"][ids]).max()
                        ),
                    }
                )
        fig, axes = plt.subplots(3, 2, figsize=(13, 8), sharex=True)
        for source in (0, 1):
            old, new = ga["predictions"][source], gb["predictions"][source]
            for i, seg in enumerate(seq["segments"]):
                j = int(seg["endIndex"])
                ids = slice(max(0, j - 60), min(n, j + 160))
                t = seq["t"][ids] - seq["t"][j - 1]
                ax = axes[i, source]
                ax.plot(
                    t,
                    old["backTorque"][ids],
                    color=AMBER,
                    lw=1.8,
                    label="Previous MPC: back",
                )
                ax.plot(
                    t,
                    new["backTorque"][ids],
                    color=BLUE,
                    lw=1.6,
                    label="Recovery MPC: back",
                )
                ax.plot(
                    t,
                    new["frontTorque"][ids],
                    color="#8054a6",
                    lw=1,
                    label="Recovery MPC: abdomen",
                )
                ax.plot(
                    t,
                    new["assistanceTorque"][ids],
                    color=GREEN,
                    ls=":",
                    lw=1.3,
                    label="Recovery assistance",
                )
                ax.axvline(0, color="#777777", lw=0.8)
                ax.axhline(12, color="#555555", ls="--", lw=0.8)
                ax.set_title(
                    f"{seg['loadKg']} kg / speed code {seg['speedCode']} | "
                    + ("Experimental" if source == 0 else "MAP"),
                    fontsize=10,
                )
                ax.set_ylabel("Torque (N·m)")
                ax.grid(alpha=0.2)
        for ax in axes[-1]:
            ax.set_xlabel("Time from recorded action end (s)")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9)
        fig.suptitle(f"Recovery transition: all endings in sequence {gb['group']}")
        fig.tight_layout(rect=(0, 0.06, 1, 0.96))
        figures.append((fig, f"mpc_recovery_endings_group{gb['group']}"))
    fig, axes = plt.subplots(2, 2, figsize=(13, 7.5))
    for source in (1, 2):
        selected = [r for r in rows if r["source"] == source]
        x = np.arange(len(selected))
        ax = axes[0, source - 1]
        ax.bar(
            x - 0.18,
            [r["beforeBackPeakNm"] for r in selected],
            0.36,
            color=AMBER,
            label="Previous MPC",
        )
        ax.bar(
            x + 0.18,
            [r["afterBackPeakNm"] for r in selected],
            0.36,
            color=BLUE,
            label="Recovery MPC",
        )
        ax.set_ylabel("Post-action back peak (N·m)")
        ax.set_title("Experimental activation" if source == 1 else "MAP activation")
        ax = axes[1, source - 1]
        ax.bar(
            x - 0.18,
            [r["beforeRecoveryAngleRMSE"] for r in selected],
            0.36,
            color=AMBER,
        )
        ax.bar(
            x + 0.18, [r["afterRecoveryAngleRMSE"] for r in selected], 0.36, color=BLUE
        )
        ax.set_ylabel("First-second standing RMSE (rad)")
        for ax in axes[:, source - 1]:
            ax.set_xticks(
                x, [f"{r['loadKg']}kg\n{r['speedCode']}" for r in selected], fontsize=8
            )
            ax.grid(axis="y", alpha=0.2)
    axes[0, 0].legend(frameon=False)
    fig.suptitle(
        "One controller for all nine conditions | unchanged Hill parameters and references"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    figures.insert(0, (fig, "mpc_recovery_summary"))
    names = []
    with PdfPages(out / "mpc_recovery_before_after.pdf") as pdf:
        for fig, name in figures:
            for suffix in ("png", "pdf"):
                filename = name + "." + suffix
                fig.savefig(out / filename, dpi=170, bbox_inches="tight")
                names.append(filename)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    names.append("mpc_recovery_before_after.pdf")
    with (out / "mpc_recovery_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    names.append("mpc_recovery_metrics.csv")
    summary = {
        "beforeWorstPeakNm": max(r["beforeBackPeakNm"] for r in rows),
        "afterWorstPeakNm": max(r["afterBackPeakNm"] for r in rows),
        "beforeMeanPeakNm": float(np.mean([r["beforeBackPeakNm"] for r in rows])),
        "afterMeanPeakNm": float(np.mean([r["afterBackPeakNm"] for r in rows])),
        "window": "First second after each recorded action end; windows truncated at sequence end",
        "records": rows,
    }
    (out / "mpc_recovery_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    names.append("mpc_recovery_summary.json")
    target.mkdir(exist_ok=True)
    for name in names:
        shutil.copy2(out / name, target / name)
    print(json.dumps({k: v for k, v in summary.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
