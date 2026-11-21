"""Prepare complete archived 5/10/15-kg sequences, including original idle seams.

No reference rescale or new waiting values are applied. The hand-written
Combination_raw already contains theta=1.5 and torque=20 during sewn waits;
these are filled reference values, not independent standing measurements.
"""
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.io import loadmat, savemat
from scipy.signal import lfilter

from prepare_hill_trials import ROOT, GROUPS, CONDITIONS, CONTACT, read_sto
from package_hill_results import struct_row


def prepare(raw_root):
    manifest = json.loads((ROOT / "config/workflows/training_map__muscle_tendon__general_model_op_scale.json").read_text())
    records, audit = [], []
    global_offset = 0
    for group_id, group in enumerate(GROUPS):
        raw_path = ROOT / f"data/processed/training_map__muscle_tendon__{group}_done__build_data/Combination_raw.mat"
        raw = loadmat(raw_path)["Combination"]
        entry = next(item for item in manifest["files"] if item["target"] == f"Combination_{group}.mat")
        mapped = loadmat(ROOT / entry["source"])["Combination"]
        assert np.array_equal(raw[:, 2], mapped[:, 2])
        assert np.array_equal(raw[:, 0], mapped[:, 3])
        samples = np.arange(0, len(raw), 10)
        n = len(samples)
        front = np.zeros(n)
        nominal = np.zeros(n)
        speed_codes = np.zeros(n)
        stage_id = np.zeros(n)
        attached = np.zeros(n)
        segments = []
        offset = 0
        waiting = np.ones(n)
        for stage, (load, speed) in enumerate(CONDITIONS[group_id*3:group_id*3+3]):
            base = raw_root / f"dynamics_stoop_{load}kg_{speed}"
            lmt, activation = [read_sto(base / name) for name in ("lmt.sto", "Activations.sto")]
            start, end = (49 if stage == 2 else 0), (349 if stage == 1 else 299)
            count = int(np.floor((lmt[end, 0]-lmt[start, 0])/.001+.0003))+1
            indices = np.flatnonzero((samples >= offset) & (samples < offset+count))
            timestamps = lmt[start, 0]+(samples[indices]-offset)*.001
            front[indices] = np.interp(timestamps, activation[:, 0], activation[:, 1])
            nominal[indices] = load
            speed_codes[indices] = speed
            stage_id[indices] = stage+1
            active = indices[(raw[samples[indices], 1] > 1e-10) | (mapped[samples[indices], 1] > 1e-10)]
            assert len(active) > 10
            waiting[active[0]:active[-1]+1] = 0
            left, right = CONTACT[group_id*3+stage]
            attached[indices] = load*((samples[indices]+global_offset+1 > left) &
                                      (samples[indices]+global_offset+1 < right))
            segments.append({"loadKg": load, "speedCode": speed,
                             "startIndex": active[0]+1, "endIndex": active[-1]+1,
                             "startTime": samples[active[0]]*.001, "endTime": samples[active[-1]]*.001})
            offset += count
        assert offset == len(raw)
        decay = np.exp(-.01/.08)
        mass = lfilter(1-decay, [1, -decay], attached, zi=[decay*attached[0]])[0]
        decay = np.exp(-.01/.12)
        front = lfilter(1-decay, [1, -decay], front, zi=[decay*front[0]])[0]
        q, torque = raw[samples, 2].copy(), raw[samples, 0].copy()
        neutral = float(np.median(q[samples < 1000]))
        assert np.isclose(neutral, 1.5, atol=1e-12, rtol=0)
        sewn = (raw[samples, 1] == 0) & np.isclose(q, neutral, atol=1e-12, rtol=0) & np.isclose(torque, 20, atol=1e-12, rtol=0)
        records.append({
            "name": group, "t": samples*.001, "q": q, "torque": torque,
            "a": np.c_[raw[samples, 1], mapped[samples, 1]], "front": front,
            "m": mass, "loadKg": nominal, "speedCode": speed_codes, "stageId": stage_id,
            "waiting": waiting, "referenceIsFilled": sewn.astype(float),
            "standingAngle": neutral, "v0": 0.0, "dt": .01,
            "segments": struct_row(segments),
            "referencePolicy": "Exact original full Combination_raw at 100 Hz; original wait seams unchanged",
        })
        audit.append({"group": group, "samples": n, "durationSeconds": samples[-1]*.001,
                      "standingAngleRad": neutral, "referenceMaximumDifference": 0.,
                      "waitTorqueIsOriginalFilledValue": 20.})
        global_offset += len(raw)
    cache = ROOT / "data/processed/hill_compliant/paired_sequences.mat"
    savemat(cache, {"sequences": struct_row(records),
                   "referenceAudit": struct_row(audit)}, do_compression=True)
    print("Saved continuous sequences:", cache)
    for row in audit:
        print(row)
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path,
                        default=ROOT.parent / "YOUHAN LI code and project description/Data/240528Pilot")
    prepare(parser.parse_args().raw_root)
