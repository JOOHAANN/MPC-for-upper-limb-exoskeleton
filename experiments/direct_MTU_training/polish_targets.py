"""Additional joint torque/angle search for routes still below fitting goals."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.optimize import least_squares
import fit_direct_MTU as fit


def task(label):
    name,s=label.split('_');source=int(s)-1;start=time.time()
    previous=json.loads((fit.HERE/f'refined_{label}.json').read_text())
    original=json.loads((fit.HERE/f'fit_{label}.json').read_text())
    (fit.HERE/f'pre_polish_{label}.json').write_text(json.dumps(previous,indent=2))
    subject=loadmat(fit.HERE/'before_training/data'/f'subject_{name}.mat',simplify_cells=True)['subject']
    fit.Q_POINT=1.5;fit.Q_SHAPE=2.;fit.T_POINT=1.5;fit.T_SHAPE=3.
    fit.RATIO_LIMIT=.25;fit.RATIO_WEIGHT=4.;fit.PEAK_RATIO_WEIGHT=.5
    fit.Q_TARGET=.705 if source==0 else .605;fit.T_TARGET=.805 if source==0 else .755
    fit.Q_GOAL_WEIGHT=40.;fit.T_GOAL_WEIGHT=60.
    x0=fit.encode(previous['parameters']);x1=fit.encode(original['parameters'])
    rng=np.random.default_rng(int(name)+source)
    rr=fit.records(subject,source,.02)
    candidates=[x0,x1,.5*x0+.5*x1,.75*x0+.25*x1]
    candidates += [np.clip(x0+rng.normal(0,.025,28),1e-6,1-1e-6) for _ in range(8)]
    def merit(x):
        residual,_=fit.evaluate(fit.decode(x),rr);return float(residual@residual)
    candidates.sort(key=merit);x=candidates[0];calls=0
    for dt,iterations in [(.02,25),(.01,15)]:
        rr=fit.records(subject,source,dt);best=[np.inf,x.copy()]
        def objective(z):
            nonlocal calls
            calls+=1;residual,_=fit.evaluate(fit.decode(z),rr);score=float(residual@residual)
            if score<best[0]:best[:]=[score,z.copy()]
            return residual
        objective(x)
        result=least_squares(objective,x,bounds=(np.zeros(28),np.ones(28)),max_nfev=iterations,
                             diff_step=1e-4,ftol=5e-5,xtol=5e-5,gtol=1e-4,x_scale='jac')
        score,x=best;_,metrics=fit.evaluate(fit.decode(x),rr)
        record=dict(subject=name,source=source+1,dt=dt,objective=score,metrics=metrics,
                    parameters=fit.decode(x).tolist(),evaluations=calls,seconds=time.time()-start,
                    status=int(result.status),message=str(result.message),
                    objectiveWeights=dict(thetaPoint=1.5,thetaShape=2,torquePoint=1.5,torqueShape=3,
                        abdomenMeanRatioSoftPreference=.25,abdomenRatioWeight=4,peakRatioWeight=.5,
                        meanThetaTarget=fit.Q_TARGET,meanTorqueTarget=fit.T_TARGET,
                        thetaGoalWeight=40,torqueGoalWeight=60))
        (fit.HERE/f'polished_{label}_{dt:.2f}.json').write_text(json.dumps(record,indent=2))
        print(label,dt,score,metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
    (fit.HERE/f'polished_{label}.json').write_text(json.dumps(record,indent=2))
    (fit.HERE/f'refined_{label}.json').write_text(json.dumps(record,indent=2))
    savemat(fit.HERE/f'refined_{label}.mat',dict(parameters=fit.decode(x),source=source+1,subject=name))
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--routes',nargs='+',default=['725_1','730_1','730_2']);args=p.parse_args()
    with ProcessPoolExecutor(max_workers=len(args.routes)) as pool:
        for future in as_completed([pool.submit(task,label) for label in args.routes]):
            r=future.result();print('DONE',r['subject'],r['source'],r['metrics'],flush=True)
