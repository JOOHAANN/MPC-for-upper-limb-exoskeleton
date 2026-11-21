"""Refine force tracking and non-dominant abdomen under the same fixed physics."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.optimize import least_squares
import fit_direct_MTU as fit


def task(job):
    name,source=job;label=f'{name}_{source+1}';start=time.time()
    initial=json.loads((fit.HERE/f'fit_{label}.json').read_text())
    subject=loadmat(fit.HERE/'before_training/data'/f'subject_{name}.mat',simplify_cells=True)['subject']
    fit.Q_POINT=1.5;fit.Q_SHAPE=2.;fit.T_POINT=1.5;fit.T_SHAPE=3.
    fit.RATIO_LIMIT=.25;fit.RATIO_WEIGHT=4.;fit.PEAK_RATIO_WEIGHT=.5
    fit.Q_TARGET=.7 if source==0 else .6;fit.T_TARGET=.8 if source==0 else .75
    x=fit.encode(initial['parameters']);calls=0
    for dt,iterations in [(.02,30),(.01,18)]:
        rr=fit.records(subject,source,dt);best=[np.inf,x.copy()]
        def objective(z):
            nonlocal calls
            calls+=1;residual,_=fit.evaluate(fit.decode(z),rr);score=float(residual@residual)
            if score<best[0]:best[:]=[score,z.copy()]
            return residual
        objective(x)
        result=least_squares(objective,x,bounds=(np.zeros(28),np.ones(28)),max_nfev=iterations,
                             diff_step=1e-4,ftol=8e-5,xtol=8e-5,gtol=1e-4,x_scale='jac')
        score,x=best;_,metrics=fit.evaluate(fit.decode(x),rr)
        record=dict(subject=name,source=source+1,dt=dt,objective=score,metrics=metrics,
                    parameters=fit.decode(x).tolist(),evaluations=calls,seconds=time.time()-start,
                    status=int(result.status),message=str(result.message),
                    objectiveWeights=dict(thetaPoint=1.5,thetaShape=2,torquePoint=1.5,torqueShape=3,
                        abdomenMeanRatioSoftPreference=.25,abdomenRatioWeight=4,peakRatioWeight=.5,
                        meanThetaTarget=fit.Q_TARGET,meanTorqueTarget=fit.T_TARGET))
        (fit.HERE/f'refined_{label}_{dt:.2f}.json').write_text(json.dumps(record,indent=2))
        print(label,dt,score,metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
    savemat(fit.HERE/f'refined_{label}.mat',dict(parameters=fit.decode(x),source=source+1,subject=name))
    (fit.HERE/f'refined_{label}.json').write_text(json.dumps(record,indent=2))
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--subjects',nargs='+',default=['528','617','725','730'])
    p.add_argument('--workers',type=int,default=8);args=p.parse_args()
    jobs=[(name,source) for name in args.subjects for source in range(2)]
    missing=[f'{name}_{source+1}' for name,source in jobs if not(fit.HERE/f'fit_{name}_{source+1}.json').exists()]
    if missing:raise RuntimeError('Initial fits not yet available: '+', '.join(missing))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for future in as_completed([pool.submit(task,job) for job in jobs]):
            r=future.result();print('DONE',r['subject'],r['source'],r['metrics'],flush=True)
