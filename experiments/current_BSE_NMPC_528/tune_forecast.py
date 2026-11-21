"""Select a causal activation-residual predictor by leave-one-condition-out CV."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.io import loadmat,savemat
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
file=ROOT/'Demo/BSE/Direct_Activation_No_Limits/data/subject_528.mat'
subject=loadmat(file,simplify_cells=True)['subject'];points=40
def matrices(lags):
    xx=[];yy=[];folds=[];ww=[]
    for group,q in enumerate(subject['sequences']):
        observed=np.c_[q['backActivation'][:,0],q['abdomenActivation']]
        planned=np.c_[q['backActivation'][:,1],q['abdomenMAPActivation']]
        residual=observed-planned;residual[q['waiting']>.5]=0;n=len(residual)
        for stage,seg in enumerate(q['segments']):
            current=np.arange(max(0,int(seg['startIndex'])-41),int(seg['endIndex']),2)
            X=np.column_stack([residual[np.maximum(0,current-lag),channel] for lag in lags for channel in [0,1]])
            Y=np.column_stack([residual[np.minimum(n-1,current+k),channel] for k in range(1,points+1) for channel in [0,1]])
            slope=(observed[current,0]-observed[np.maximum(0,current-4),0])/.04
            weights=1+4*np.clip(slope/.5,0,1)
            xx.append(X);yy.append(Y);ww.append(weights);folds.extend([3*group+stage]*len(current))
    return np.vstack(xx),np.vstack(yy),np.array(folds),np.concatenate(ww)

reports=[];best=None
for lags in [np.arange(0,41,2),np.arange(0,41),np.r_[np.arange(0,21),25,30,40]]:
    X,Y,folds,w=matrices(lags)
    for penalty in [1e-6,1e-5,1e-4,.001,.01,.1]:
        errors=[]
        for fold in range(9):
            train=folds!=fold;test=~train;scale=np.maximum(X[train].std(0),.002)
            Z=X[train]/scale;W=w[train]
            beta=np.linalg.solve(Z.T@(W[:,None]*Z)+penalty*np.eye(Z.shape[1]),Z.T@(W[:,None]*Y[train]))
            error=X[test]/scale@beta-Y[test]
            weights=np.tile(np.repeat(np.linspace(2,1,points),2),(len(error),1))
            errors.append(float(np.sqrt(np.average(np.mean(error**2*weights,axis=1),weights=w[test]))))
        score=float(np.mean(errors));reports.append(dict(lags=lags.tolist(),penalty=penalty,weightedCVRMSE=score,foldRMSE=errors))
        if best is None or score<best[0]:best=(score,lags.copy(),penalty)
score,lags,penalty=best;X,Y,folds,w=matrices(lags);scale=np.maximum(X.std(0),.002);Z=X/scale
beta=np.linalg.solve(Z.T@(w[:,None]*Z)+penalty*np.eye(Z.shape[1]),Z.T@(w[:,None]*Y))
predictor=dict(lags=lags.astype(float),scale=scale,coefficients=beta,points=40,predictionDt=.01,penalty=penalty,
    usesOnlyPastActivationResiduals=True,usesAngleOrTorqueLabels=False,
    calibrationRMSE=float(np.sqrt(np.mean((Z@beta-Y)**2))),baselineMAPResidualRMSE=float(np.sqrt(np.mean(Y**2))),
    weightedCVRMSE=score,calibrationPolicy='Leave-one-condition-out lag/ridge selection; extra weight on rising activation samples')
modelSHA256=hashlib.sha256(file.read_bytes()).hexdigest()
savemat(ROOT/'Demo/NMPC/data/current_activation_predictor_528.mat',dict(predictor=predictor,modelSHA256=modelSHA256),do_compression=True,long_field_names=True)
(HERE/'activation_predictor_CV.json').write_text(json.dumps(dict(selected=dict(lags=lags.tolist(),penalty=penalty,weightedCVRMSE=score),reports=reports),indent=2))
print(json.dumps(dict(lags=lags.tolist(),penalty=penalty,weightedCVRMSE=score,calibrationRMSE=predictor['calibrationRMSE']),indent=2))
