"""Select completed strict candidates using all action samples, not plot shifts."""
import csv,json
import numpy as np
import strict_fit_528 as f

subject=f.loadmat(f.HERE/'strict_input_528.mat',simplify_cells=True)['subject']
selected=[];table=[]
for source in range(2):
    rr=f.records(subject,source,.01);candidates=[]
    paths=[f.HERE/f'strict_528_{source+1}.json',f.HERE/f'strict_global_528_{source+1}.json',f.HERE/f'strict_polished_528_{source+1}.json']
    for path in paths:
        if not path.exists():continue
        report=json.loads(path.read_text())
        if report['dt']!=.01:continue
        p=np.ascontiguousarray(report['parameters']);gain=report['abdomenScale'];qr=[];tr=[];qe=[];te=[];ratios=[];rests=[]
        complete=True;bounds=[];frontPeaks=[];backPeaks=[]
        for q,r in zip(f.cells(subject['sequences']),rr):
            out=f.rollout(p,gain,r);complete &= bool(np.all(np.isfinite(out)))
            moving=r['waiting']<.5
            if complete:
                assert np.max(np.abs(out[moving,19:22]))==0
                assert np.max(np.abs(out[:,23:26]))==0
            rests.append(float(out[0,4]))
            for seg in f.cells(q['segments']):
                ix=np.arange(int(seg['startIndex'])-1,int(seg['endIndex']));valid=ix[q['torqueValid'][ix].astype(bool)]
                a=out[ix,0];b=q['referenceTheta'][ix];t=out[valid,2];ref=q['referenceTorque'][valid]
                qr.append(f.corr(a,b));tr.append(f.corr(t,ref))
                qe.append(float(np.rad2deg(np.sqrt(np.mean((a-b)**2)))));te.append(float(np.sqrt(np.mean((t-ref)**2))))
                ratios.append(float(np.mean(out[ix,4])/max(np.mean(out[ix,3]),1)))
                bounds.extend([float(a.min()),float(a.max())]);frontPeaks.append(float(out[ix,4].max()));backPeaks.append(float(out[ix,3].max()))
        targetT=.8 if source==0 else .75;targetQ=.7 if source==0 else .6
        metrics=dict(ThetaR=float(np.mean(qr)),TorqueR=float(np.mean(tr)),ThetaRMSEDeg=float(np.mean(qe)),
            TorqueRMSENm=float(np.mean(te)),AbdomenBackMeanRatio=float(np.mean(ratios)),StandingAbdomenNm=float(np.mean(rests)),
            MinimumThetaRad=min(bounds),MaximumThetaRad=max(bounds),AbdomenPeakNm=max(frontPeaks),BackPeakNm=max(backPeaks))
        passed=complete and metrics['ThetaR']>targetQ and metrics['TorqueR']>targetT
        # Prefer small abdominal recruitment once fitting targets are met.
        score=metrics['ThetaRMSEDeg']/15+metrics['TorqueRMSENm']/40+8*max(metrics['AbdomenBackMeanRatio']-.2,0)+max(metrics['StandingAbdomenNm']-3,0)/3
        row=dict(Source=source+1,File=path.name,Complete=complete,TargetsMet=passed,SelectionScore=score,**metrics)
        table.append(row);candidates.append((not passed,score,path.name))
    candidates.sort();assert candidates and not candidates[0][0],f'No completed source {source+1} fit meets both mean-correlation targets.'
    selected.append(candidates[0][2])
selection=dict(subject='528',files=selected,usesAllActionSamples=True,
    policy='First satisfy mean torque/theta correlation targets. Then minimize normalized RMSE plus soft abdominal/recruitment preferences.',
    score='theta_RMSE_deg/15 + torque_RMSE_Nm/40 + 8*max(mean_abdomen_back_ratio-0.2,0) + max(initial_standing_abdomen_Nm-3,0)/3',
    independentValidation=False)
(f.HERE/'strict_selection.json').write_text(json.dumps(selection,indent=2))
with (f.HERE/'strict_candidate_comparison.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
print(json.dumps(selection,indent=2));print(json.dumps(table,indent=2))
