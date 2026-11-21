"""Additional multistart search for the strict abdominal-input calibration."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,threading,argparse
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.optimize import least_squares
import strict_fit_528 as f

def main(source):
    start=time.time();subject=f.loadmat(f.HERE/'strict_input_528.mat',simplify_cells=True)['subject']
    rr=f.records(subject,source,.04);rng=np.random.default_rng(528+source);seeds=[]
    # Sample broad effective anatomical/activation initializations. Reference
    # MT lengths are formed from positive slack lengths, never angle limits.
    for k in range(240):
        p=f.cells(subject['models'])[source]['parameters'].copy()
        p[0]=rng.uniform(150,1500);p[8]=rng.uniform(10,240)
        p[1]=rng.uniform(.25,.9);p[2]=rng.uniform(.65,1.1)
        p[3]=rng.uniform(5,45);p[4]=rng.uniform(3,12);p[5]=rng.uniform(1.2,3)
        p[6]=rng.uniform(.025,.09);p[7]=rng.uniform(.015,.18)
        p[9]=rng.uniform(.4,1.8);p[10]=rng.uniform(.25,.8);p[11]=rng.uniform(.55,.95)
        p[12]=rng.uniform(.35,.8);p[13]=rng.uniform(.1,.9)
        p[14]=rng.uniform(20,45);p[15]=rng.uniform(.2,.35)
        p[16]=rng.uniform(.025,.06);p[17]=rng.uniform(.035,.1)
        p[20]=rng.uniform(3,12);p[21]=rng.uniform(.025,.09);p[22]=rng.uniform(.015,.18);p[23]=rng.uniform(.1,.3)
        p[18]=rng.uniform(.05,.3)+p[16]/p[1]*np.sqrt(p[2]**2-np.sin(p[23])**2)
        p[19]=rng.uniform(.05,.35)+p[17]/p[10]*np.sqrt(p[11]**2-np.sin(p[23])**2)
        seeds.append(f.encode(p,rng.uniform(.1,4)))
    def evaluate(z):
        p,g=f.decode(z);r,m=f.evaluate(p,g,rr,source);return float(r@r),m
    lock=threading.Lock();calls=0
    with ThreadPoolExecutor(max_workers=8) as pool:
        scored=list(pool.map(evaluate,seeds));order=np.argsort([s[0] for s in scored]);bestOverall=[np.inf,seeds[order[0]]]
        print('SOURCE',source+1,'GLOBAL SEEDS',[(i,round(scored[i][0],2),scored[i][1]['thetaR'],scored[i][1]['torqueR'],scored[i][1]['abdomenBackRatio']) for i in order[:5]],flush=True)
        for seedIndex in range(3):
            rr=f.records(subject,source,.04);x=seeds[order[seedIndex]];best=[np.inf,x.copy()]
            def objective(z):
                nonlocal calls
                p,g=f.decode(z);r,m=f.evaluate(p,g,rr,source);value=float(r@r)
                with lock:
                    calls+=1
                    if value<best[0]:best[:]=[value,z.copy()]
                return r
            result=least_squares(objective,x,bounds=(np.zeros(25),np.ones(25)),max_nfev=55,
                diff_step=1e-4,ftol=1e-4,xtol=1e-4,gtol=1e-4,x_scale='jac',workers=pool.map)
            p,g=f.decode(best[1]);_,metrics=f.evaluate(p,g,rr,source)
            print('SOURCE',source+1,'START',seedIndex,'SCORE',best[0],metrics,'seconds',round(time.time()-start),flush=True)
            if best[0]<bestOverall[0]:bestOverall=best.copy()
            rec=dict(source=source+1,dt=.04,parameters=p.tolist(),abdomenScale=g,objective=best[0],metrics=metrics)
            (f.HERE/f'strict_global_528_{source+1}_seed{seedIndex}.json').write_text(json.dumps(rec,indent=2))
        for dt in [.02,.01]:
            rr=f.records(subject,source,dt);x=bestOverall[1];best=[np.inf,x.copy()]
            result=least_squares(objective,x,bounds=(np.zeros(25),np.ones(25)),max_nfev=60,
                diff_step=8e-5,ftol=4e-5,xtol=4e-5,gtol=1e-4,x_scale='jac',workers=pool.map)
            p,g=f.decode(best[1]);_,metrics=f.evaluate(p,g,rr,source);bestOverall=best.copy()
            rec=dict(subject='528',source=source+1,dt=dt,parameters=p.tolist(),abdomenScale=g,objective=best[0],metrics=metrics,
                     inputPolicy='Recorded abdomen or its independently trained MAP; one amplitude gain only',
                     movementPD=False,standingPD=True,angleBoundaries=False,independentValidation=False,
                     evaluations=calls,seconds=time.time()-start)
            (f.HERE/f'strict_global_528_{source+1}.json').write_text(json.dumps(rec,indent=2))
            print('SOURCE',source+1,'REFINE',dt,best[0],metrics,'seconds',round(time.time()-start),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('source',type=int);args=parser.parse_args();main(args.source-1)
