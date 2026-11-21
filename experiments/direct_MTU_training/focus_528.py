"""Prioritize subject 528 with parallel numerical derivatives and multistart."""
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor,as_completed
import json
from pathlib import Path
import threading
import time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.optimize import least_squares
import fit_direct_MTU as fit


def task(source):
    label=f'528_{source+1}';start=time.time()
    current=json.loads((fit.HERE/f'refined_{label}.json').read_text())
    original=json.loads((fit.HERE/f'fit_{label}.json').read_text())
    (fit.HERE/f'pre_focus_{label}.json').write_text(json.dumps(current,indent=2))
    subject=loadmat(fit.HERE/'before_training/data/subject_528.mat',simplify_cells=True)['subject']
    fit.Q_POINT=1.8;fit.Q_SHAPE=3.;fit.T_POINT=1.;fit.T_SHAPE=2.
    fit.RATIO_LIMIT=.25;fit.RATIO_WEIGHT=4.;fit.PEAK_RATIO_WEIGHT=.5
    fit.Q_TARGET=.80;fit.T_TARGET=.85;fit.Q_GOAL_WEIGHT=40.;fit.T_GOAL_WEIGHT=60.
    rr=fit.records(subject,source,.02)
    def residual_for(z):
        p=fit.decode(z);residual,metrics=fit.evaluate(p,rr);extra=[]
        # Fit the two bending depths and excessive overshoot in the objective
        # only. These terms never apply forces or constrain simulated states.
        for r in rr:
            out=fit.rollout(p,r)
            for ix,_ in r['cases']:
                pred=np.nan_to_num(out[ix,0],nan=10);ref=r['theta'][ix];mid=len(ix)//2
                extra += [(.7/.2)*(np.min(pred[:mid])-np.min(ref[:mid])),
                          (.7/.2)*(np.min(pred[mid:])-np.min(ref[mid:])),
                          (.7/.2)*max(np.max(pred)-np.max(ref),0)]
        return np.r_[residual,extra],metrics
    x0=fit.encode(current['parameters']);x1=fit.encode(original['parameters'])
    seeds=[x0,x1,.5*x0+.5*x1,.75*x0+.25*x1]
    rng=np.random.default_rng(5280+source)
    seeds += [np.clip(x0+rng.normal(0,.015,28),1e-6,1-1e-6) for _ in range(4)]
    def score(x):r,_=residual_for(x);return float(r@r)
    seeds.sort(key=score);best_overall=[np.inf,seeds[0].copy()];calls=0
    lock=threading.Lock()
    with ThreadPoolExecutor(max_workers=8) as pool:
        for seed_index,x in enumerate(seeds[:2]):
            stage_best=[np.inf,x.copy()]
            def objective(z):
                nonlocal calls
                residual,_=residual_for(z);value=float(residual@residual)
                with lock:
                    calls+=1
                    if value<stage_best[0]:stage_best[:]=[value,z.copy()]
                return residual
            objective(x)
            result=least_squares(objective,x,bounds=(np.zeros(28),np.ones(28)),
                max_nfev=30,diff_step=1e-4,ftol=5e-5,xtol=5e-5,gtol=1e-4,
                x_scale='jac',workers=pool.map)
            value,x=stage_best;_,metrics=residual_for(x)
            if value<best_overall[0]:best_overall[:]=[value,x.copy()]
            print(label,'seed',seed_index,value,metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
        rr=fit.records(subject,source,.01);x=best_overall[1];stage_best=[np.inf,x.copy()]
        objective(x)
        result=least_squares(objective,x,bounds=(np.zeros(28),np.ones(28)),
            max_nfev=22,diff_step=8e-5,ftol=3e-5,xtol=3e-5,gtol=1e-4,
            x_scale='jac',workers=pool.map)
        value,x=stage_best;_,metrics=residual_for(x)
    record=dict(subject='528',source=source+1,dt=.01,objective=value,metrics=metrics,
                parameters=fit.decode(x).tolist(),evaluations=calls,seconds=time.time()-start,
                status=int(result.status),message=str(result.message),
                objectiveWeights=dict(thetaPoint=1.8,thetaShape=3,torquePoint=1,torqueShape=2,
                    abdomenMeanRatioSoftPreference=.25,abdomenRatioWeight=4,peakRatioWeight=.5,
                    meanThetaTarget=.8,meanTorqueTarget=.85,thetaGoalWeight=40,torqueGoalWeight=60,
                    extremumWeight=.7,extremumScaleRad=.2))
    (fit.HERE/f'focused_{label}.json').write_text(json.dumps(record,indent=2))
    (fit.HERE/f'refined_{label}.json').write_text(json.dumps(record,indent=2))
    savemat(fit.HERE/f'refined_{label}.mat',dict(parameters=fit.decode(x),source=source+1,subject='528'))
    print(label,'FINAL',value,metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
    return record


if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(task,s) for s in range(2)]):future.result()
