"""Calibrate causal activation-residual forecasts; no torque/angle labels."""
from pathlib import Path
import sys, json
import numpy as np
from scipy.io import loadmat, savemat

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
sys.path.insert(0,str(ROOT/'experiments/controller_comparison_528'))
from tune_pid import inputs

def main():
    subject=loadmat(ROOT/'Demo/BSE/data/subject_528.mat',simplify_cells=True)['subject']
    summaries=[]
    for name,lags in [('history_1s',np.r_[np.arange(0,41,2),np.arange(45,101,5)]),
                      ('history_3s',np.r_[np.arange(0,41,2),np.arange(45,101,5),np.arange(110,301,10)])]:
        xx=[]; yy=[]; gg=[]
        for group,seq in enumerate(subject['sequences']):
            model=subject['models'][0]
            observed=inputs(seq,model,0); planned=inputs(seq,model,1)
            residual=observed-planned; count=len(observed)
            selected=np.zeros(count,dtype=bool)
            for seg in seq['segments']:
                selected[max(0,int(seg['startIndex'])-31):min(count-42,int(seg['endIndex'])+50)]=True
            current=np.flatnonzero(selected)[::2]
            history=np.maximum(0,current[:,None]-lags[None,:])
            x=residual[history].reshape(len(current),-1)
            query=current[:,None]+np.arange(1,17)[None,:]*.025/seq['dt']
            y=np.stack([np.interp(query,np.arange(count),residual[:,j]) for j in (0,1)],axis=-1).reshape(len(current),32)
            xx.append(x); yy.append(y); gg.extend([group]*len(x))
        x=np.vstack(xx); y=np.vstack(yy); groups=np.array(gg)
        scale=np.maximum(np.sqrt(np.mean(x*x,axis=0)),1e-4); z=x/scale
        trials=[]
        for penalty in [1e-5,1e-4,.001,.01,.1,1]:
            scores=[]
            for fold in range(3):
                train=groups!=fold; test=~train; a=z[train]
                c=np.linalg.solve(a.T@a+penalty*np.sum(train)*np.eye(z.shape[1]),a.T@y[train])
                scores.append(float(np.sqrt(np.mean((z[test]@c-y[test])**2))))
            trials.append(dict(penalty=penalty,leaveOneSequenceOutRMSE=float(np.mean(scores)),foldRMSE=scores))
        best=min(trials,key=lambda r:r['leaveOneSequenceOutRMSE']); penalty=best['penalty']
        c=np.linalg.solve(z.T@z+penalty*len(z)*np.eye(z.shape[1]),z.T@y)
        predictor=dict(lags=lags.astype(float),scale=scale,coefficients=c,points=16.,predictionDt=.025,
            penalty=penalty,usesOnlyPastActivationResiduals=True,usesAngleOrTorqueLabels=False,
            calibrationRMSE=float(np.sqrt(np.mean((z@c-y)**2))),
            baselineMAPResidualRMSE=float(np.sqrt(np.mean(y*y))),
            leaveOneSequenceOutRMSE=best['leaveOneSequenceOutRMSE'])
        savemat(HERE/f'forecast_{name}.mat',dict(predictor=predictor),long_field_names=True)
        result=dict(name=name,featureCount=z.shape[1],samples=len(z),**{k:v for k,v in predictor.items() if np.isscalar(v)},crossValidation=trials)
        summaries.append(result); print(json.dumps(result),flush=True)
    (HERE/'forecast_training.json').write_text(json.dumps(summaries,indent=2)+'\n')

if __name__=='__main__': main()
