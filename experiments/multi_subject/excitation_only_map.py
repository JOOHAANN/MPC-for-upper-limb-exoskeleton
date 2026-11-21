"""Build MAP labels and complete bouts from excitation alone: no subject IK.

Only emgFilt.sto is read from each recording. A=0.1 supplies activation labels;
the frozen first-subject smoothing envelope is a generic signal benchmark.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from scipy.io import savemat

from audit_subjects import DATA, ROOT, read_sto
from quality_gate import gaussian_matlab, sha
from review_exp_one_cycle import find_m_bouts

CHANNELS = [
    p + s for p in ("iliocost", "longlumb", "longthor") for s in ("_r", "_l")
] + ["rectabd_r", "rectabd_l"]


def excitation_signal(path):
    """Read normalized excitation and its native clock, without any geometry."""
    names, table = read_sto(path)
    if not np.isfinite(table).all() or np.any(np.diff(table[:, 0]) < 0):
        raise ValueError("Invalid excitation clock or values")
    values = table[:, [names.index(n) for n in CHANNELS]]
    times, first, counts = np.unique(table[:, 0], return_index=True, return_counts=True)
    conflicts = []
    for i in np.flatnonzero(counts > 1):
        repeated = values[table[:, 0] == times[i]]
        if not np.allclose(repeated, repeated[0], atol=1e-10, rtol=1e-9):
            conflicts.append(float(times[i]))
    values = values[first]
    allowed = ~np.isin(times, conflicts)
    times, values = times[allowed], values[allowed]
    fine = np.arange(0, times[-1] - times[0], 0.001)
    sampled = np.column_stack([np.interp(fine, times - times[0], v) for v in values.T])
    # Do not silently rectify a signed export. Range checks remain local to
    # the selected bout, and unusable bouts remain explicit exclusions.
    excitation = gaussian_matlab(gaussian_matlab(sampled, 5), 300)[::10]
    activation = gaussian_matlab(
        gaussian_matlab(np.expm1(0.1 * sampled) / np.expm1(0.1), 5), 300
    )[::10]
    return dict(
        timeEpochSeconds=times[0] + fine[::10],
        excitation=excitation,
        activation=activation,
        nativeClock=times,
        excludedIntervals=[(v - 0.2, v + 0.2) for v in conflicts],
        conflictingTimes=conflicts,
    )


def build_review(subject, output):
    output = Path(output) / subject
    output.mkdir(parents=True, exist_ok=True)
    frozen = json.loads(
        (ROOT / "docs/subject_quality_v2/quality_gate.json").read_text()
    )
    limits = frozen["smoothnessUpperLimits"]
    choices = {}
    audits = []
    for path in sorted(
        (DATA / subject / "CalibrationData").glob("dynamics_stoop*/emgFilt.sto")
    ):
        match = re.fullmatch(
            r"dynamics_stoop_(\d+(?:\.\d+)?)kg_(\d+)(?:_(\d+))?", path.parent.name
        )
        if not match:
            continue
        load, speed = float(match[1]), int(match[2])
        if load not in (5, 7.5, 10, 12.5, 15):
            continue
        trial = str(path.parent.relative_to(DATA / subject))
        try:
            s = excitation_signal(path)
            t = s["timeEpochSeconds"]
            six = s["activation"][:, :6]
            bouts = find_m_bouts(
                t,
                six.mean(1),
                six,
                [s["nativeClock"]],
                limits,
                excluded_intervals=s["excludedIntervals"],
            )
            valid = []
            for b in bouts:
                ids = (t >= b["startEpochSeconds"]) & (t <= b["endEpochSeconds"])
                if (
                    s["activation"][ids].min() < -1e-6
                    or s["activation"][ids].max() > 1 + 1e-6
                ):
                    b["concerns"].append("activation_outside_zero_one")
                if not b["concerns"] and not b["smoothnessWarnings"]:
                    valid.append(b)
            audit = dict(
                trial=trial,
                loadKg=load,
                speedCode=speed,
                sourceSHA256=sha(path),
                completeBouts=len(bouts),
                acceptedBouts=len(valid),
                usesIK=False,
                usesCEINMSTorque=False,
            )
            if valid:
                b = valid[0]
                key = (load, speed)
                rank = (bool(match[3]), b["startEpochSeconds"])
                if key not in choices or rank < choices[key][0]:
                    ids = (t >= b["startEpochSeconds"]) & (t <= b["endEpochSeconds"])
                    target = output / f"{load:g}kg_{speed}.mat"
                    record = dict(
                        subject=subject.removesuffix("_A01"),
                        trial=trial,
                        loadKg=load,
                        speedCode=speed,
                        selected=b,
                        windowPath=str(target),
                    )
                    choices[key] = (
                        rank,
                        record,
                        {
                            k: s[k][ids]
                            for k in ("timeEpochSeconds", "excitation", "activation")
                        },
                    )
            audits.append(audit)
        except Exception as exc:
            audits.append(
                dict(trial=trial, loadKg=load, speedCode=speed, error=str(exc))
            )
    conditions = []
    for key, (_, record, window) in sorted(choices.items()):
        savemat(record["windowPath"], window, do_compression=True)
        conditions.append(record)
    payload = dict(
        version="excitation-only-map-v1",
        subject=subject,
        conditions=conditions,
        trials=audits,
        shapeFactor=0.1,
        requiresIK=False,
        requiresMuscleForces=False,
        requiresCEINMSTorque=False,
        smoothnessUpperLimits=limits,
    )
    audit_file = output / "excitation_review.json"
    audit_file.write_text(json.dumps(payload, indent=2) + "\n")
    gate = dict(
        version=payload["version"],
        smoothnessUpperLimits=limits,
        excitationAuditPath=str(audit_file),
        excitationAuditSHA256=sha(audit_file),
    )
    return gate, payload
