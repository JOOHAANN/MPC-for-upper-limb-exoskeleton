"""Export exact metrics, fixed load scaling and input provenance for subject 528."""
import csv,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.io import loadmat
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
DEMO=ROOT/'Demo/BSE/Direct_Activation_No_Limits'
subject=loadmat(DEMO/'data/subject_528.mat',simplify_cells=True)['subject']
mapData=loadmat(DEMO/'data/legacy_MAP_528.mat',simplify_cells=True)
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
rows=list(csv.DictReader((DEMO/'output/case_metrics.csv').open()));assert len(rows)==18
summary=[]
for label in ['Experimental activation','MAP activation']:
    selected=[r for r in rows if r['Source']==label];assert len(selected)==9
    avg=lambda key:float(np.mean([float(r[key]) for r in selected]))
    peak=lambda key:max(float(r[key]) for r in selected)
    assert all(r['Complete'].lower() in ['1','true'] for r in selected)
    record=dict(Source=label,Cases=9,TorqueCorrelation=avg('TorqueCorrelation'),ThetaCorrelation=avg('ThetaCorrelation'),
        TorqueRMSENm=avg('TorqueRMSENm'),ThetaRMSEDeg=avg('ThetaRMSEDeg'),AbdomenBackMeanRatio=avg('AbdomenBackMeanRatio'),
        AbdomenPeakNm=peak('AbdomenPeakNm'),BackPeakNm=peak('BackPeakNm'),
        MinimumThetaRad=min(float(r['MinimumThetaRad']) for r in selected),MaximumThetaRad=peak('MaximumThetaRad'),
        MovementPDPeakNm=peak('MovementPDPeakNm'),BoundaryPeakNm=peak('BoundaryPeakNm'))
    assert record['MovementPDPeakNm']==0 and record['BoundaryPeakNm']==0
    assert record['TorqueCorrelation']>(.8 if label.startswith('Experimental') else .75)
    assert record['ThetaCorrelation']>(.7 if label.startswith('Experimental') else .6)
    summary.append(record)
def write_csv(path,records):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
write_csv(DEMO/'output/summary.csv',summary)
s=subject['MAPScaling'];scaling=[]
for i,weight in enumerate(s['loadKg']):
    for j,muscle in enumerate(['Back','Abdomen']):
        scaling.append(dict(LoadKg=float(weight),Muscle=muscle,Gain=float(s['gain'][i,j]),Offset=float(s['offset'][i,j]),
            ExperimentalMeanMinimum=float(s['experimentalMinimum'][i,j]),ExperimentalMeanMaximum=float(s['experimentalMaximum'][i,j]),
            RawMAPMeanMinimum=float(s['rawMAPMinimum'][i,j]),RawMAPMeanMaximum=float(s['rawMAPMaximum'][i,j]),
            SpeedCodes='25/35/45/55/65'))
write_csv(DEMO/'data/MAP_load_scaling.csv',scaling)

rawProvenance=[];selection=[]
for r in mapData['records']:
    path=Path(r['sourcePath'])
    for name in ['Activations.sto','emgFilt.sto','lmt.sto','MuscleForces.sto','ma_L5_S1_Flex_Ext.sto']:
        rawProvenance.append(dict(path=str(path/name),sha256=sha(path/name)))
    selection.append(dict(LoadKg=float(r['loadKg']),SpeedCode=int(r['speedCode']),
        StartEpoch=float(r['startEpoch']),EndEpoch=float(r['endEpoch']),Selection=r['selection'],SourceFolder=str(path)))
write_csv(DEMO/'data/experimental_windows.csv',selection)
previous=json.loads((ROOT/'experiments/controller_comparison_528/original_file_hashes.json').read_text())
protected={k:v for k,v in previous.items() if k.startswith('Demo\\BSE\\')}
assert len(protected)==20 and all(sha(ROOT/k)==v for k,v in protected.items())
oldAudit=json.loads((ROOT/'experiments/direct_MTU_training/abdominal_MAP_528/training_audit.json').read_text())
assert all(sha(Path(r['sourcePath']))==r['sourceSHA256'] for r in oldAudit['provenance'])
legacy=Path(subject['MAP']['backSource']);legacyFiles=['Primitives_fit_y.mat','baseline_data.mat','pv_coefficients.mat','pw_coefficients.mat',
    'MAP_model.m','getWeights.m','MAP_output_model_fusion.m','Dynamic_kine_reference.m','data_combination.m','NNMF_component.m']
manifest=dict(createdUTC=datetime.now(timezone.utc).isoformat(),subject='528',modelVersion=subject['modelVersion'],
    backMAPPolicy='Reuse original trained predictor and preserve the repeated IL term literally',
    abdominalMAPPolicy=subject['MAP']['abdomen']['architecture'],physicalCoefficientsPerSource=24,
    physicalInputGainsOptimized=0,loadScalingFrozenBeforePhysicalCalibration=True,scaling=scaling,
    experimentalSpeedCodes=[25,35,45,55,65],referenceCurvesModified=False,
    movementPD=False,standingPD=True,assistanceNm=0,angleBoundaries=False,independentBSEValidation=False,
    rawFiles=rawProvenance,legacyFiles={name:sha(legacy/name) for name in legacyFiles},
    verifiedOriginalBSEFilesUnchanged=20,verifiedOriginalAbdominalSourceFilesUnchanged=25,
    inputPackageSHA256=sha(HERE/'calibration_input.mat'),publishedPackageSHA256=sha(DEMO/'data/subject_528.mat'),
    outputSHA256=sha(DEMO/'output/BSE_results.mat'),summary=summary,pendingRevisedSubjects=['617','725','730'])
(DEMO/'data/source_manifest.json').write_text(json.dumps(manifest,indent=2))
report=['# Subject 528: legacy selection and load-only MAP rescaling','',
    'Original back MAP coefficients and the repeated IL fusion term are retained. Abdominal MAP was independently retrained.',
    'For each load and MTU, five normalized-phase experimental profiles are averaged first; their global minimum/maximum set the target amplitude range.',
    'One fixed affine mapping per load/MTU applies to every speed. No extra input gain is optimized during physical calibration.',
    'Two independent 24-parameter physical vectors were recalibrated over the same nine bouts, one per source.',
    'Movement PD, assistance and boundary torque are zero. Standing recruitment remains enabled outside recorded actions.','',
    '| Source | Mean torque r | Mean theta r | Mean abdomen/back torque ratio |',
    '|---|---:|---:|---:|']
for r in summary:report.append(f"| {r['Source']} | {r['TorqueCorrelation']:.4f} | {r['ThetaCorrelation']:.4f} | {r['AbdomenBackMeanRatio']:.3f} |")
report+=['','Metrics are arithmetic means of nine action-only coefficients, using all recorded action samples.',
    'These are fitting replays, not held-out BSE validation. Small abdomen/back ratios are model outputs, not independently measured physiological ratios.',
    'The original trial-specific 50:150-point central-valley correction is replaced by the shared load-only global-range mapping.',
    'The 15-kg/code-35 angle fit remains the weakest condition; see the complete per-condition table.',
    'The old 5-kg/code-25 and code-35 EMG windows exceed CEINMS export support. Earlier complete M pairs are retained from the same recordings; no missing activation is extrapolated.',
    'All original reference curves, timelines, mass inputs and standing settings were verified unchanged by the MATLAB release test.',
    'The revised method currently applies to 528 only. Subjects 617/725/730 remain pending.','',
    'Read summary.csv, case_metrics.csv, data/MAP_load_scaling.csv and data/experimental_windows.csv for exact values.']
(DEMO/'output/528_method_and_results.md').write_text('\n'.join(report)+'\n')
print(json.dumps(summary,indent=2));print('PASS: mean targets, zero movement PD/boundaries, source hashes, legacy-source manifest and fixed load-scaling export.')
