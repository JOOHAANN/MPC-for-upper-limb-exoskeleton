"""Write a provenance manifest and a compact report for the strict 528 demo."""
import csv,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import strict_fit_528 as f

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
rows=list(csv.DictReader((f.DEMO/'output/case_metrics.csv').open()))
assert len(rows)==18 and all(r['Complete']=='1' or r['Complete'].lower()=='true' for r in rows)
selection=json.loads((f.HERE/'strict_selection.json').read_text())
summary=[]
for source,label in enumerate(['Experimental activation','MAP activation']):
    rr=[r for r in rows if r['Source']==label];assert len(rr)==9
    avg=lambda key:float(np.mean([float(r[key]) for r in rr]))
    peak=lambda key:max(float(r[key]) for r in rr)
    report=json.loads((f.HERE/selection['files'][source]).read_text())
    row=dict(Subject='528',Source=label,Cases=9,TorqueCorrelation=avg('TorqueCorrelation'),
        ThetaCorrelation=avg('ThetaCorrelation'),TorqueRMSENm=avg('TorqueRMSENm'),ThetaRMSEDeg=avg('ThetaRMSEDeg'),
        AbdomenBackMeanRatio=avg('AbdomenBackMeanRatio'),AbdomenPeakNm=peak('AbdomenPeakNm'),
        BackPeakNm=peak('BackPeakNm'),AbdomenScale=report['abdomenScale'],
        MovementPDPeakNm=peak('MovementPDPeakNm'),BoundaryPeakNm=peak('BoundaryPeakNm'))
    assert row['TorqueCorrelation']>(.8 if source==0 else .75)
    assert row['ThetaCorrelation']>(.7 if source==0 else .6)
    assert row['MovementPDPeakNm']==0 and row['BoundaryPeakNm']==0
    summary.append(row)
with (f.DEMO/'output/summary.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)

previous=json.loads((f.HERE.parent/'controller_comparison_528/original_file_hashes.json').read_text())
original={k:v for k,v in previous.items() if k.startswith('Demo\\BSE\\')}
assert len(original)==20 and all(sha(f.HERE.parents[1]/k)==v for k,v in original.items())
mapAudit=json.loads((f.HERE/'abdominal_MAP_528/training_audit.json').read_text())
assert len(mapAudit['provenance'])==25
for record in mapAudit['provenance']:
    assert sha(Path(record['sourcePath']))==record['sourceSHA256']
files=[
    'data/subject_528.mat','data/abdominal_MAP_528.mat','prepare_activation.m',
    'predict_abdominal_MAP.m','dual_MTU_model.m','private/native/direct_mtu_mex.c',
    'private/native/mtu_mechanics.h','output/BSE_results.mat','output/case_metrics.csv']
manifest=dict(createdUTC=datetime.now(timezone.utc).isoformat(),strictSubjects=['528'],
    pendingStrictSubjects=['617','725','730'],selectedTrainingReports=selection['files'],
    recordedAbdominalChannel='rect_abd_r',experimentalPreprocessing='Original stitching, 100-Hz interpolation, causal 0.12-s filter',
    abdominalMAPArchitecture=mapAudit['architecture'],MAPUsesIK=False,MAPUsesTorque=False,MAPUsesBSEOutput=False,
    parameterPolicy='One 24-parameter physical vector plus one abdominal amplitude gain per source; nine bouts jointly fitted',
    movementPD=False,standingPD=True,angleBoundaries=False,assistanceNm=0,
    outerStepSeconds=.01,internalStepSeconds=.001,independentBSEValidation=False,
    originalBSEFilesVerifiedUnchanged=20,originalAbdominalFilesVerifiedUnchanged=25,
    hashes={name:sha(f.DEMO/name) for name in files},summary=summary)
(f.DEMO/'data/source_manifest.json').write_text(json.dumps(manifest,indent=2))
report=['# Strict abdominal-input BSE: subject 528','',
    'All nine conditions use one fixed physical vector and abdominal amplitude gain per source.',
    'Experimental abdomen uses recorded rect_abd_r only. MAP abdomen is trained separately from that experimental channel.',
    'Standing feedback is retained; movement feedback, assistance and angle boundaries are zero.',
    'References and back-MAP outputs are unchanged. Metrics below average the nine action-only Pearson coefficients.',
    'These are calibration replays, not held-out BSE validation. Effective fitted parameters are not independently identified anatomical measurements.','',
    '| Input | Torque r | Theta r | Abdomen/back mean torque ratio | Abdomen gain |',
    '|---|---:|---:|---:|---:|']
for r in summary:
    report.append(f"| {r['Source']} | {r['TorqueCorrelation']:.4f} | {r['ThetaCorrelation']:.4f} | {r['AbdomenBackMeanRatio']:.3f} | {r['AbdomenScale']:.4f} |")
report+=['','The abdominal/back ratio is the arithmetic mean of per-action mean-torque ratios, not a measured physiological ratio.',
    'Weakest angle fit: 15 kg / code 35, experimental r = 0.4633 and MAP r = 0.4462.',
    'MAP minimum angle is approximately -0.095 rad; it is retained without clipping.',
    'Abdominal MAP held-out-load correlation before final refitting is approximately 0.33; its small activation RMSE does not establish strong shape generalization.',
    'Old mixed-abdomen candidates are archived outside the active demo. Subjects 617, 725 and 730 still require this strict-input refit.','',
    'See summary.csv, case_metrics.csv and data/trained_parameters.csv for exact values.']
(f.DEMO/'output/strict_528_report.md').write_text('\n'.join(report)+'\n')
print(json.dumps(manifest,indent=2))
