"""Publish the ablation honestly, including controller and model limitations."""
from pathlib import Path
import hashlib,json,csv
from datetime import datetime,timezone
import numpy as np
from scipy.io import loadmat
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent.parent
OUT=ROOT/'Demo/NMPC/output/no_human_PD_528'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    prior=json.loads((HERE/'prior_NMPC_hashes.json').read_text())
    protected=json.loads((ROOT/'experiments/controller_comparison_528/original_file_hashes.json').read_text())
    frozen=json.loads((ROOT/'results/nine_case_75pct/release_manifest.json').read_text())
    frozen={**frozen['artifactSHA256'],**frozen['sourceCodeSHA256']}
    for p,h in {**prior,**protected,**frozen}.items():assert digest(ROOT/p)==h,p
    log=(HERE/'verification_final.log').read_text(errors='replace')
    assert 'PASS: no human movement/standing PD or integral torque' in log
    r=loadmat(OUT/'comparison_no_human_PD_528.mat',simplify_cells=True)['results']
    with (OUT/'case_metrics.csv').open(newline='',encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    assert len(rows)==54
    summary=[]
    for source in ('Experimental activation','MAP activation'):
        for method in ('Unassisted','PID','NMPC'):
            selected=[row for row in rows if row['Source']==source and row['Controller']==method]
            def values(key):return np.array([float(row[key]) for row in selected])
            summary.append(dict(Source=source,Controller=method,Conditions=len(selected),
                MeanExcessPeakNm=float(np.mean(values('ExcessPeakNm'))),
                TotalExcessAreaNms=float(np.sum(values('ExcessAreaNms'))),
                NoExceedanceCases=int(np.sum(values('ExcessPeakNm')<1e-8)),
                MeanThetaTaskRMSEDeg=float(np.mean(values('ThetaTaskRMSEDeg'))),
                WorstThetaTaskErrorDeg=float(np.max(values('ThetaTaskMaximumErrorDeg'))),
                MaximumAbdomenTorqueNm=float(np.max(values('AbdomenPeakNm'))),
                MaximumBoundaryTorqueMagnitudeNm=float(np.max(values('SoftBoundaryPeakMagnitudeNm'))),
                MeanAbsolutePositiveDemandLagMs=float(np.nanmean(np.abs(values('PositiveDemandLagMs'))))
                    if np.any(np.isfinite(values('PositiveDemandLagMs'))) else None))
    with (OUT/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    contract=[]
    for group,g in enumerate(r['groups'],1):
        for source,pair in enumerate(g['runs'],1):
            for method in ('pid','nmpc'):
                actual=pair[method]['actual']
                flags=np.atleast_1d(pair[method]['solverExitFlag'])
                times=np.atleast_1d(pair[method]['solveSeconds']); times=times[times>0]
                contract.append(dict(group=group,source=source,controller=method,
                    maximumHumanPostureTerm=float(np.max(np.abs(actual[:,21:23]))),
                    minimumTorqueNm=float(np.min(actual[:,24])),maximumTorqueNm=float(np.max(actual[:,24])),
                    minimumThetaRad=float(np.min(actual[:,0])),maximumThetaRad=float(np.max(actual[:,0])),
                    failedOptimizationUpdates=int(np.sum(np.isfinite(flags)&(flags<=0))),
                    meanControllerComputeSeconds=float(np.mean(times)),
                    p95ControllerComputeSeconds=float(np.percentile(times,95))))
    report=['# Subject-528 no-human-PD torque/trajectory control','',
        'The original comparison and calibrated model files are unchanged. Movement posture Kp/Kd and standing Kp/Kd/Ki are zero in the actual and predictive human plants. All sampled posture correction and integral outputs are zero. PID and NMPC both act only through the same signed external actuator.','',
        'Both controllers receive the same fixed MAP task trajectory. It is an offline unassisted rollout of the previously calibrated MAP model: original posture gains exist only in this reference-generation step. They do not run inside the ablated actual or prediction plant. No future experimental activation or IK is passed to either controller. This known-task preview must not be described as online unknown-intention estimation.','',
        'Both controllers use signed assistance in [-80, +80] N m, torque-magnitude growth at 200 N m/s, immediate withdrawal, and 20-ms updates. Reversal withdraws to zero before opposite torque builds. The previous one-way hardware assumption has therefore changed, with explicit user authorization. Squared torque effort is not matched.','',
        'The PID includes a reactive muscle-torque branch and an external angle PID, calibrated across the recorded cases. NMPC predicts 20 points at 0.025 s with four blocks over 0.5 s and jointly penalizes task-angle/velocity errors, clipped MAP-task torque tracking and threshold exceedance. Targets compete when actuation or the retained human excitation response limits control authority. Thresholds retain 75% of each successful-task unassisted peak from the preceding calibration, equally for both controllers. The drifting no-PD ablation is not used to redefine the task threshold.','',
        '| Source | Controller | Mean peak excess (N m) | Total excess area (N m s) | Cases within threshold | Mean task-angle RMSE (deg) | Max abdomen torque (N m) | Max boundary torque (N m) |',
        '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for s in summary:
        report.append(f"| {s['Source']} | {s['Controller']} | {s['MeanExcessPeakNm']:.2f} | {s['TotalExcessAreaNms']:.2f} | {s['NoExceedanceCases']}/9 | {s['MeanThetaTaskRMSEDeg']:.2f} | {s['MaximumAbdomenTorqueNm']:.2f} | {s['MaximumBoundaryTorqueMagnitudeNm']:.2f} |")
    report+=['','## Interpretation','',
        'Angle error is evaluated against the declared task trajectory, not against the drifting uncontrolled no-PD motion. A separate CSV column records deviation from that uncontrolled motion. These two errors answer different questions. All RQ1/RQ2 figures use the same raw closed-loop trajectories, with no output clipping, smoothing or per-case gain changes.','',
        'Removing the posture gains changes the calibrated human model structurally. Its previous fit and validation do not carry over. Gravity compensation, passive damping and the original soft boundaries remain. The retained assistance-dependent excitation response also remains: during standing, its requested muscle torque contains -gravity - assistance. It can therefore substantially cancel sustained external assistance. This limits posture-control authority even with signed torque and helps explain imperfect tracking; no replacement hidden human PD is added.','',
        'Positive-demand lag is reported separately because active braking and trajectory control can require signed assistance when torque is below its threshold. Smaller angle error alone does not prove improved temporal alignment or reduced physical hardware latency. There is no equal-effort comparison and no hardware/comfort validation.','',
        'All nine cases per source are calibration/evaluation on the recorded protocol, not independent validation. Controller/model contract checks pass; success on a research question must be judged from the performance table, not a QP success flag.','',
        f'Preservation audit: {len(prior)} prior NMPC runtime files, {len(protected)} original Demo files, and {len(frozen)} frozen-release files match their previous hashes.','',
        'Run `Demo/NMPC/offline_PID_NMPC_no_human_PD_528.m`. Results stay in `Demo/NMPC/output/no_human_PD_528/`. Standard all-nine torque/theta, continuous examples, and abdominal/soft-boundary torque audit figures are provided in PNG/PDF/FIG.']
    (OUT/'comparison_report.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    result=dict(createdAtUTC=datetime.now(timezone.utc).isoformat(),contractChecksPassed=True,
        all18CasesEvaluated=True,previousNMPCFilesUnchanged=len(prior),originalDemoFilesUnchanged=len(protected),
        frozenReleaseFilesUnchanged=len(frozen),humanPosturePDEnabled=False,
        matchedActuatorLimits=[-80,80],magnitudeRiseRateNmPerSec=200,effortMatched=False,
        independentValidation=False,referenceProvenance=r['protocol']['thetaReference'],
        contract=contract,summary=summary)
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
