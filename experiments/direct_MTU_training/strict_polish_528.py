"""Prioritize the requested mean correlations without changing runtime physics."""
import argparse,json,time,threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.optimize import least_squares
import strict_fit_528 as f

def main(source):
    start=time.time();subject=f.loadmat(f.HERE/'strict_input_528.mat',simplify_cells=True)['subject']
    rr=f.records(subject,source,.02);targetQ=.725 if source==0 else .67;targetT=.83 if source==0 else .79
    def loss(z):
        p,g=f.decode(z);r,m=f.evaluate(p,g,rr,source)
        r[-2]=120*max(targetQ-m['thetaR'],0)
        r[-1]=80*max(targetT-m['torqueR'],0)
        return r,m
    files=list(f.HERE.glob(f'strict_global_528_{source+1}_seed*.json'))
    files+=list(f.HERE.glob(f'strict_global_528_{source+1}.json'))
    files+=list(f.HERE.glob(f'strict_528_{source+1}.json'))
    candidates=[]
    for path in files:
        record=json.loads(path.read_text());z=f.encode(record['parameters'],record['abdomenScale']);r,m=loss(z)
        candidates.append((float(r@r),z,path.name))
    candidates.sort(key=lambda c:c[0]);x=candidates[0][1]
    print('SOURCE',source+1,'POLISH FROM',candidates[0][2],flush=True)
    lock=threading.Lock();calls=0
    with ThreadPoolExecutor(max_workers=8) as pool:
        for dt in [.02,.01]:
            rr=f.records(subject,source,dt);best=[np.inf,x.copy()]
            def objective(z):
                nonlocal calls
                r,m=loss(z);value=float(r@r)
                with lock:
                    calls+=1
                    if value<best[0]:best[:]=[value,z.copy()]
                return r
            objective(x)
            result=least_squares(objective,x,bounds=(np.zeros(25),np.ones(25)),max_nfev=65,
                diff_step=6e-5,ftol=2e-5,xtol=2e-5,gtol=1e-4,x_scale='jac',workers=pool.map)
            x=best[1];p,g=f.decode(x);_,metrics=loss(x)
            record=dict(subject='528',source=source+1,dt=dt,parameters=p.tolist(),abdomenScale=g,
                objective=best[0],metrics=metrics,evaluations=calls,seconds=time.time()-start,
                inputPolicy='Recorded abdomen or independently trained abdominal MAP, constant gain only',
                movementPD=False,standingPD=True,angleBoundaries=False,independentValidation=False,
                thetaGoal=targetQ,thetaGoalWeight=120,torqueGoal=targetT,torqueGoalWeight=80)
            (f.HERE/f'strict_polished_528_{source+1}.json').write_text(json.dumps(record,indent=2))
            print('SOURCE',source+1,'POLISHED',dt,best[0],metrics,'seconds',round(time.time()-start),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source',type=int);args=parser.parse_args();main(args.source-1)
