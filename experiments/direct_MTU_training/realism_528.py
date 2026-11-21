"""528 refinement with explicit low-resting-coactivation and recovery checks.

These are calibration penalties, not torque caps, angle boundaries or runtime
feedback. The forward simulator remains unchanged and receives no references.
"""
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor,as_completed
import argparse
import json
import threading
import time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.optimize import least_squares
import fit_direct_MTU as fit


def task(source):
    label=f'528_{source+1}';start=time.time()
    initial=json.loads((fit.HERE/f'refined_{label}.json').read_text())
    (fit.HERE/f'pre_realism_{label}.json').write_text(json.dumps(initial,indent=2))
    subject=loadmat(fit.HERE/'before_training/data/subject_528.mat',simplify_cells=True)['subject']
    fit.Q_POINT=1.5;fit.Q_SHAPE=2.5;fit.T_POINT=1.2;fit.T_SHAPE=2.5
    fit.RATIO_LIMIT=.25;fit.RATIO_WEIGHT=6.;fit.PEAK_RATIO_WEIGHT=.5
    fit.Q_TARGET=.75;fit.T_TARGET=.85;fit.Q_GOAL_WEIGHT=40.;fit.T_GOAL_WEIGHT=60.
    rr=fit.records(subject,source,.02)
    def loss(z):
        p=fit.decode(z);outputs=[fit.rollout(p,r) for r in rr]
        residual,metrics=fit.evaluate(p,rr,outputs=outputs);extra=[];rest=[]
        for r,out in zip(rr,outputs):
            rest.append(float(out[0,4]))
            # The standing controller requests 2 N m coactivation. A broad
            # 5-N-m tonic abdominal preference discourages passive-force
            # cancellation at rest; it is not an experimentally measured cap.
            extra.append(max(np.nan_to_num(out[0,4],nan=1000)-5,0)/(5*np.sqrt(len(rr))))
            for ix,_ in r['cases']:
                pred=np.nan_to_num(out[ix,0],nan=10);ref=r['theta'][ix];mid=len(ix)//2
                extra += [(.5/.2)*(np.min(pred[:mid])-np.min(ref[:mid])),
                          (.5/.2)*(np.min(pred[mid:])-np.min(ref[mid:])),
                          (.5/.2)*max(np.max(pred)-np.max(ref),0)]
                post=np.arange(ix[-1]+1,min(len(out),ix[-1]+round(1.2/r['dt'])+1))
                post=post[r['waiting'][post]>.5]
                if len(post):
                    q=np.nan_to_num(out[post,0],nan=10)
                    v=np.nan_to_num(out[post,1],nan=10)
                    extra.extend((q-r['cfg'][0])/(.2*np.sqrt(len(post))))
                    extra.extend(.2*v/np.sqrt(len(post)))
        metrics['standingAbdomenTorqueNm']=float(np.mean(rest))
        return np.r_[residual,extra],metrics
    x0=fit.encode(initial['parameters']);seeds=[x0]
    for back_length,front_length in [(None,.9),(.9,.9),(1.,.9),(1.,.95)]:
        x=x0.copy()
        if back_length is not None:x[2]=(back_length-fit.LOW[2])/(fit.HIGH[2]-fit.LOW[2])
        x[11]=(front_length-fit.LOW[11])/(fit.HIGH[11]-fit.LOW[11]);seeds.append(x)
    def score(x):r,_=loss(x);return float(r@r)
    seeds.sort(key=score);best_total=[np.inf,seeds[0].copy()];calls=0;lock=threading.Lock()
    with ThreadPoolExecutor(max_workers=8) as pool:
        for seed_index,x in enumerate(seeds[:2]):
            stage_best=[np.inf,x.copy()]
            def objective(z):
                nonlocal calls
                residual,_=loss(z);value=float(residual@residual)
                with lock:
                    calls+=1
                    if value<stage_best[0]:stage_best[:]=[value,z.copy()]
                return residual
            objective(x)
            result=least_squares(objective,x,bounds=(np.zeros(28),np.ones(28)),
                max_nfev=35,diff_step=1e-4,ftol=5e-5,xtol=5e-5,gtol=1e-4,
                x_scale='jac',workers=pool.map)
            value,x=stage_best;_,metrics=loss(x)
            if value<best_total[0]:best_total[:]=[value,x.copy()]
            print(label,'seed',seed_index,value,metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
        rr=fit.records(subject,source,.01);x=best_total[1];stage_best=[np.inf,x.copy()]
        objective(x)
        result=least_squares(objective,x,bounds=(np.zeros(28),np.ones(28)),
            max_nfev=25,diff_step=8e-5,ftol=3e-5,xtol=3e-5,gtol=1e-4,
            x_scale='jac',workers=pool.map)
        value,x=stage_best;_,metrics=loss(x)
    record=dict(subject='528',source=source+1,dt=.01,objective=value,metrics=metrics,
                parameters=fit.decode(x).tolist(),evaluations=calls,seconds=time.time()-start,
                status=int(result.status),message=str(result.message),
                objectiveWeights=dict(thetaPoint=1.5,thetaShape=2.5,torquePoint=1.2,torqueShape=2.5,
                    abdomenMeanRatioSoftPreference=.25,abdomenRatioWeight=6,peakRatioWeight=.5,
                    meanThetaTarget=.75,meanTorqueTarget=.85,thetaGoalWeight=40,torqueGoalWeight=60,
                    extremumWeight=.5,extremumScaleRad=.2,restingAbdomenPreferenceNm=5,
                    restingAbdomenScaleNm=5,standingRecoverySeconds=1.2,recoveryAngleScaleRad=.2))
    (fit.HERE/f'realistic_{label}.json').write_text(json.dumps(record,indent=2))
    (fit.HERE/f'refined_{label}.json').write_text(json.dumps(record,indent=2))
    savemat(fit.HERE/f'refined_{label}.mat',dict(parameters=fit.decode(x),source=source+1,subject='528'))
    print(label,'FINAL',value,metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--sources',nargs='+',type=int,default=[1,2]);args=parser.parse_args()
    with ProcessPoolExecutor(max_workers=len(args.sources)) as pool:
        for future in as_completed([pool.submit(task,s-1) for s in args.sources]):future.result()
