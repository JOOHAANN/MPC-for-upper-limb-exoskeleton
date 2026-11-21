"""Select the completed fits using every action sample at the final 1-ms step."""
import sys,json,csv
from pathlib import Path
import numpy as np
from scipy.io import loadmat
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'direct_MTU_training'))
import strict_fit_528 as f
subject=loadmat(HERE/'calibration_input.mat',simplify_cells=True)['subject']
selected=[];comparison=[]
for source in range(2):
    rr=f.records(subject,source,.01);choices=[]
    files=[HERE/f'fit_528_{source+1}.json',HERE/f'fit_528_{source+1}_cross.json']
    for path in files:
        if not path.exists():continue
        report=json.loads(path.read_text())
        if report['dt']!=.01:continue
        p=np.ascontiguousarray(report['parameters']);qr=[];tr=[];qe=[];te=[];ratios=[];rests=[]
        complete=True
        for q,r in zip(f.cells(subject['sequences']),rr):
            out=f.rollout(p,1,r);complete &= bool(np.all(np.isfinite(out)));rests.append(float(out[0,4]))
            assert np.max(np.abs(out[:,23:26]))==0
            assert np.max(np.abs(out[r['waiting']<.5,19:22]))==0
            for seg in f.cells(q['segments']):
                ix=np.arange(int(seg['startIndex'])-1,int(seg['endIndex']));valid=ix[q['torqueValid'][ix].astype(bool)]
                qr.append(f.corr(out[ix,0],q['referenceTheta'][ix]));tr.append(f.corr(out[valid,2],q['referenceTorque'][valid]))
                qe.append(float(np.rad2deg(np.sqrt(np.mean((out[ix,0]-q['referenceTheta'][ix])**2)))))
                te.append(float(np.sqrt(np.mean((out[valid,2]-q['referenceTorque'][valid])**2))))
                ratios.append(float(out[ix,4].mean()/max(out[ix,3].mean(),1)))
        metrics=dict(ThetaR=float(np.mean(qr)),TorqueR=float(np.mean(tr)),ThetaRMSEDeg=float(np.mean(qe)),
            TorqueRMSENm=float(np.mean(te)),AbdomenBackRatio=float(np.mean(ratios)),StandingAbdomenNm=float(np.mean(rests)))
        passed=complete and metrics['ThetaR']>(.7 if source==0 else .6) and metrics['TorqueR']>(.8 if source==0 else .75)
        score=metrics['ThetaRMSEDeg']/15+metrics['TorqueRMSENm']/40+8*max(metrics['AbdomenBackRatio']-.2,0)+max(metrics['StandingAbdomenNm']-3,0)/3
        comparison.append(dict(Source=source+1,File=path.name,Complete=complete,TargetsMet=passed,SelectionScore=score,**metrics))
        choices.append((not passed,score,path.name))
    choices.sort();assert choices and not choices[0][0],'No completed fit meets the requested source-wise mean targets.'
    selected.append(choices[0][2])
result=dict(files=selected,usesAllActionSamples=True,independentValidation=False,
    policy='Satisfy source-wise mean correlation targets, then prefer lower normalized RMSE and smaller abdominal recruitment.')
(HERE/'selected_fits.json').write_text(json.dumps(result,indent=2))
with (HERE/'candidate_comparison.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(comparison[0]));writer.writeheader();writer.writerows(comparison)
print(json.dumps(result,indent=2));print(json.dumps(comparison,indent=2))
