"""Refit 24 physical coefficients per source with frozen load-only MAP scaling."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,threading,argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor,as_completed
import numpy as np
from scipy.io import loadmat
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'direct_MTU_training'))
import strict_fit_528 as f

def encode(p):return f.encode(p,1)[:24]
def decode(z):return f.decode(np.r_[z,(1-f.LOW[24])/(f.HIGH[24]-f.LOW[24])])[0]

def task(source,iterations,tag=''):
    start=time.time();f.library()
    subject=loadmat(HERE/'calibration_input.mat',simplify_cells=True)['subject']
    old=loadmat(HERE/'before_load_rescale/data/subject_528.mat',simplify_cells=True)['subject']
    rr=f.records(subject,source,.04);rng=np.random.default_rng(1528+source)
    initial=f.cells(subject['models'])[source]['parameters'];seeds=[encode(initial)]
    oldBack=np.concatenate([q['backActivation'][:,source] for q in f.cells(old['sequences'])])
    newBack=np.concatenate([q['backActivation'][:,source] for q in f.cells(subject['sequences'])])
    oldFront=np.concatenate([q['abdomenActivation' if source==0 else 'abdomenMAPActivation'] for q in f.cells(old['sequences'])])
    newFront=np.concatenate([q['abdomenActivation' if source==0 else 'abdomenMAPActivation'] for q in f.cells(subject['sequences'])])
    gain=f.cells(old['models'])[source]['abdomenScale']
    for bm,fm in [(1,1),(.75,1),(1.25,1),(1,.6),(1,1.4)]:
        p=initial.copy();p[0]*=(oldBack.mean()/max(newBack.mean(),1e-5))**p[5]*bm
        p[8]*=(gain*oldFront.mean()/max(newFront.mean(),1e-5))**p[9]*fm
        seeds.append(encode(p))
    for model in f.cells(subject['models']):
        for bm,fm in [(.8,1),(.8,1.3),(1,1.3),(1,1.6)]:
            p=model['parameters'].copy();p[0]*=bm;p[8]*=fm
            seeds.append(encode(p))
    for k in range(120):
        p=initial.copy()
        p[0]=rng.uniform(150,1500);p[8]=rng.uniform(10,400)
        p[1]=rng.uniform(.25,1.05);p[2]=rng.uniform(.65,1.1);p[3]=rng.uniform(5,50)
        p[4]=rng.uniform(3,14);p[5]=rng.uniform(.7,2.8);p[6]=rng.uniform(.025,.095);p[7]=rng.uniform(.015,.19)
        p[9]=rng.uniform(.4,1.8);p[10]=rng.uniform(.25,.9);p[11]=rng.uniform(.55,1.05)
        p[12]=rng.uniform(.35,.8);p[13]=rng.uniform(.1,.9);p[14]=rng.uniform(20,45);p[15]=rng.uniform(.2,.35)
        p[16]=rng.uniform(.025,.06);p[17]=rng.uniform(.035,.1)
        p[20]=rng.uniform(3,14);p[21]=rng.uniform(.025,.095);p[22]=rng.uniform(.015,.19);p[23]=rng.uniform(.1,.3)
        p[18]=rng.uniform(.05,.3)+p[16]/p[1]*np.sqrt(p[2]**2-np.sin(p[23])**2)
        p[19]=rng.uniform(.05,.35)+p[17]/p[10]*np.sqrt(p[11]**2-np.sin(p[23])**2)
        seeds.append(encode(p))
    label=f'528_{source+1}{tag}';path=HERE/f'fit_{label}.json'
    if path.exists():seeds.append(encode(json.loads(path.read_text())['parameters']))
    def loss(z):
        p=decode(z);r,m=f.evaluate(p,1,rr,source)
        r[-2]=70*max((.73 if source==0 else .66)-m['thetaR'],0)
        r[-1]=70*max((.83 if source==0 else .79)-m['torqueR'],0)
        return r,m
    def score(z):r,m=loss(z);return float(r@r),m
    with ThreadPoolExecutor(max_workers=8) as pool:
        scores=list(pool.map(score,seeds));order=np.argsort([x[0] for x in scores]);overall=[np.inf,seeds[order[0]]]
        print(label,'SEEDS',[(int(i),round(scores[i][0],2),scores[i][1]['thetaR'],scores[i][1]['torqueR']) for i in order[:4]],flush=True)
        lock=threading.Lock();calls=0
        def objective(z):
            nonlocal calls
            r,m=loss(z);value=float(r@r)
            with lock:
                calls+=1
                if value<best[0]:best[:]=[value,z.copy()]
            return r
        for i in range(2):
            x=seeds[order[i]];best=[np.inf,x.copy()];objective(x)
            result=least_squares(objective,x,bounds=(np.zeros(24),np.ones(24)),max_nfev=iterations,
                diff_step=1e-4,ftol=4e-5,xtol=4e-5,gtol=1e-4,x_scale='jac',workers=pool.map)
            if best[0]<overall[0]:overall=best.copy()
            _,metrics=loss(best[1]);print(label,'START',i,best[0],metrics,'seconds',round(time.time()-start),flush=True)
            rec=dict(source=source+1,dt=.04,parameters=decode(best[1]).tolist(),abdomenScale=1,metrics=metrics,objective=best[0])
            (HERE/f'fit_{label}_seed{i}.json').write_text(json.dumps(rec,indent=2))
        for dt in [.02,.01]:
            rr=f.records(subject,source,dt);x=overall[1];best=[np.inf,x.copy()];objective(x)
            result=least_squares(objective,x,bounds=(np.zeros(24),np.ones(24)),max_nfev=iterations,
                diff_step=8e-5,ftol=3e-5,xtol=3e-5,gtol=1e-4,x_scale='jac',workers=pool.map)
            overall=best.copy();_,metrics=loss(best[1])
            rec=dict(subject='528',source=source+1,dt=dt,parameters=decode(best[1]).tolist(),abdomenScale=1,
                metrics=metrics,objective=best[0],seconds=time.time()-start,evaluations=calls,
                parameterPolicy='24 physical parameters per source across nine cases; no optimized activation gains',
                scalingPolicy='One fixed affine map per load/MTU from five-speed mean experimental and raw MAP profiles',
                movementPD=False,standingPD=True,angleBoundaries=False,independentValidation=False)
            path.write_text(json.dumps(rec,indent=2));(HERE/f'fit_{label}_{dt}.json').write_text(json.dumps(rec,indent=2))
            print(label,'REFINE',dt,best[0],metrics,'seconds',round(time.time()-start),flush=True)
    return rec

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--iterations',type=int,default=60);parser.add_argument('--sources',nargs='+',type=int,default=[1,2]);parser.add_argument('--tag',default='');args=parser.parse_args()
    with ProcessPoolExecutor(max_workers=len(args.sources)) as pool:
        for future in as_completed([pool.submit(task,s-1,args.iterations,args.tag) for s in args.sources]):future.result()
