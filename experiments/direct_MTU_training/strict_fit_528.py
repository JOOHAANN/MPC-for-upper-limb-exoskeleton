"""Joint nine-bout calibration with measured/MAP abdomen and one amplitude gain.

The forward model has 24 physical parameters. Neither a reference, a decoded
movement target nor a back-derived abdominal recruitment signal enters it.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
import ctypes
import json
import threading
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor,as_completed
import numpy as np
from scipy.io import loadmat,savemat
from scipy.optimize import least_squares

HERE=Path(__file__).resolve().parent
DEMO=HERE.parents[1]/'Demo/BSE/Direct_Activation_No_Limits'
LOW=np.array([60,.15,.45,.5,2,.4,.02,.005,5,.4,.25,.45,.345,0,
              15,.20,.020,.035,.025,.025,2,.02,.005,.10,.05])
HIGH=np.array([1600,1.1,1.45,60,15,3,.10,.20,650,2,.95,1.45,.85,1.1,
               45,.35,.060,.10,.35,.40,15,.10,.20,.35,5.])
LIB=None

def library():
    global LIB
    if LIB is None:
        loaded=ctypes.CDLL(str(HERE/'direct_mtu_strict24.dll'))
        arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        loaded.direct_mtu_rollout.argtypes=[arr,arr,ctypes.c_int,arr,arr,arr,arr,ctypes.c_double,ctypes.c_int,arr]
        loaded.direct_mtu_rollout.restype=ctypes.c_int
        LIB=loaded
    return LIB

def encode(p,gain=1):
    x=np.r_[p,gain].astype(float)
    w=np.sin(x[23])
    x[18]-=x[16]/x[1]*np.sqrt(x[2]**2-w*w)
    x[19]-=x[17]/x[10]*np.sqrt(x[11]**2-w*w)
    return np.clip((x-LOW)/(HIGH-LOW),1e-5,1-1e-5)

def decode(z):
    x=LOW+z*(HIGH-LOW);p=x[:24].copy();w=np.sin(p[23])
    p[18]+=p[16]/p[1]*np.sqrt(p[2]**2-w*w)
    p[19]+=p[17]/p[10]*np.sqrt(p[11]**2-w*w)
    return np.ascontiguousarray(p),float(x[24])

def cells(x):return [x] if isinstance(x,dict) else list(x)

def records(subject,source,dt):
    rr=[]
    for q in cells(subject['sequences']):
        ix=np.arange(0,len(q['time']),round(dt/q['dt']))
        r=dict(back=q['backActivation'][ix,source],
            front=q['abdomenActivation' if source==0 else 'abdomenMAPActivation'][ix],
            mass=np.ascontiguousarray(q['mass'][ix]),waiting=np.ascontiguousarray(q['waiting'][ix],dtype=float),
            initial=np.ascontiguousarray(q['initialState'],dtype=float),
            cfg=np.ascontiguousarray(np.r_[q['standingAngle'],subject['balanceController']]),dt=dt,
            time=q['time'][ix],theta=q['referenceTheta'][ix],torque=q['referenceTorque'][ix],cases=[])
        for seg in cells(q['segments']):
            lo=q['time'][int(seg['startIndex'])-1];hi=q['time'][int(seg['endIndex'])-1]
            full=np.flatnonzero((r['time']>=lo)&(r['time']<=hi))
            chosen=full[np.unique(np.round(np.linspace(0,len(full)-1,80)).astype(int))]
            valid=chosen[q['torqueValid'][ix[chosen]].astype(bool)&np.isfinite(r['torque'][chosen])]
            post=np.arange(full[-1]+1,min(len(ix),full[-1]+round(1.2/dt)+1))
            post=post[r['waiting'][post]>.5]
            r['cases'].append((chosen,valid,post))
        rr.append(r)
    return rr

def rollout(p,gain,r):
    inp=np.ascontiguousarray(np.c_[r['back'],np.clip(gain*r['front'],0,1)])
    out=np.full((len(inp),27),np.nan)
    steps=max(1,round(r['dt']/(.002 if r['dt']>.01 else .001)))
    last=library().direct_mtu_rollout(p,r['cfg'],len(inp),inp,r['mass'],r['waiting'],r['initial'],r['dt'],steps,out)
    if last<len(inp):out[last:]=np.nan
    return out

def corr(a,b):
    if len(a)<3 or min(np.std(a),np.std(b))<1e-10:return 0.
    return float(np.corrcoef(a,b)[0,1])

def evaluate(p,gain,rr,source):
    residual=[];qr=[];tr=[];qe=[];te=[];ratios=[];rests=[];ranges=[];complete=True
    for r in rr:
        out=rollout(p,gain,r);rests.append(float(np.nan_to_num(out[0,4],nan=1000)))
        residual.append(max(rests[-1]-5,0)/5/np.sqrt(len(rr)))
        for ix,valid,post in r['cases']:
            complete &= bool(np.all(np.isfinite(out[ix,:])))
            q=np.nan_to_num(out[ix,0],nan=10,posinf=10,neginf=-10);ref=r['theta'][ix];n=np.sqrt(len(ix))
            residual.extend(1.5*(q-ref)/(.3*n))
            residual.extend(2.5*((q-q.mean())/max(q.std(),.12)-(ref-ref.mean())/max(ref.std(),.12))/n)
            qr.append(corr(q,ref));qe.append(float(np.sqrt(np.mean((q-ref)**2))));ranges.append([float(q.min()),float(q.max())])
            mid=len(ix)//2
            residual.extend([.5*(q[:mid].min()-ref[:mid].min())/.2,.5*(q[mid:].min()-ref[mid:].min())/.2,
                             .5*max(q.max()-ref.max(),0)/.2])
            t=np.nan_to_num(out[valid,2],nan=3000,posinf=3000,neginf=-3000);ref=r['torque'][valid];n=np.sqrt(len(valid))
            residual.extend(1.2*(t-ref)/(50*n))
            residual.extend(2.5*((t-t.mean())/max(t.std(),15)-(ref-ref.mean())/max(ref.std(),15))/n)
            tr.append(corr(t,ref));te.append(float(np.sqrt(np.mean((t-ref)**2))))
            back=np.nan_to_num(out[ix,3],nan=1000);front=np.nan_to_num(out[ix,4],nan=1000)
            ratio=front.mean()/max(back.mean(),1);ratios.append(float(ratio))
            residual.extend([8*max(ratio-.25,0),max(front.max()/max(back.max(),1)-.8,0)])
            if len(post):
                residual.extend((np.nan_to_num(out[post,0],nan=10)-r['cfg'][0])/(.2*np.sqrt(len(post))))
                residual.extend(.2*np.nan_to_num(out[post,1],nan=10)/np.sqrt(len(post)))
    residual.extend([.08*np.log((p[0]/p[16])/6000),.08*np.log((p[8]/p[17])/1500)])
    residual.extend([40*max((.72 if source==0 else .65)-np.mean(qr),0),60*max((.82 if source==0 else .78)-np.mean(tr),0)])
    metrics=dict(thetaR=float(np.mean(qr)),torqueR=float(np.mean(tr)),thetaRMSEDeg=float(np.rad2deg(np.mean(qe))),
        torqueRMSE=float(np.mean(te)),abdomenBackRatio=float(np.mean(ratios)),standingAbdomenTorqueNm=float(np.mean(rests)),
        minimumThetaRad=float(np.min(ranges)),maximumThetaRad=float(np.max(ranges)),complete=complete,
        caseThetaR=qr,caseTorqueR=tr)
    return np.asarray(residual),metrics

def seeds(subject,source):
    p=cells(subject['models'])[source]['parameters'];choices=[]
    for back,front,damp,gain in [(1,.2,35,1),(.3,.15,35,1),(.5,.15,45,.5),(.7,.3,45,1),(.4,.1,25,2)]:
        trial=p.copy();trial[0]*=back;trial[8]*=front;trial[3]=damp
        choices.append(encode(trial,gain))
    # Passive abdominal tension is suppressed by shorter reference fibers.
    # All values here are starting guesses, not runtime muscle or angle caps.
    for length in [.7,.85,1.]:
        for ratio in [.4,.7]:
            trial=p.copy();trial[[0,8,3,5,9,7,22]]=[400,100,35,1.,.7,.08,.08]
            trial[[2,11,10]]=[.95,length,ratio]
            trial[18]=.15+trial[16]/trial[1]*np.sqrt(trial[2]**2-np.sin(trial[23])**2)
            trial[19]=.2+trial[17]/trial[10]*np.sqrt(trial[11]**2-np.sin(trial[23])**2)
            choices.append(encode(trial,1))
    return choices

def task(source,iterations):
    start=time.time();label=f'528_{source+1}';subject=loadmat(HERE/'strict_input_528.mat',simplify_cells=True)['subject']
    rr=records(subject,source,.04)
    candidates=seeds(subject,source)
    checkpoint=HERE/f'strict_{label}.json'
    if checkpoint.exists():
        old=json.loads(checkpoint.read_text());candidates.append(encode(old['parameters'],old['abdomenScale']))
    def score(z):
        p,g=decode(z);r,m=evaluate(p,g,rr,source);return float(r@r),m
    with ThreadPoolExecutor(max_workers=8) as pool:
        values=list(pool.map(score,candidates))
        order=np.argsort([v[0] for v in values]);print(label,'SEEDS',[(round(values[i][0],2),values[i][1]['thetaR'],values[i][1]['torqueR']) for i in order[:3]],flush=True)
        best=[np.inf,candidates[order[0]]];calls=0;lock=threading.Lock()
        def objective(z):
            nonlocal calls
            p,g=decode(z);r,_=evaluate(p,g,rr,source);value=float(r@r)
            with lock:
                calls+=1
                if value<best[0]:best[:]=[value,z.copy()]
            return r
        for stage,(dt,nfev) in enumerate([(.04,iterations),(.02,iterations),(.01,iterations)]):
            rr=records(subject,source,dt);x=best[1].copy();best=[np.inf,x];objective(x)
            result=least_squares(objective,x,bounds=(np.zeros(25),np.ones(25)),max_nfev=nfev,
                diff_step=1e-4 if dt>.01 else 6e-5,ftol=4e-5,xtol=4e-5,gtol=1e-4,x_scale='jac',workers=pool.map)
            p,g=decode(best[1]);_,metrics=evaluate(p,g,rr,source)
            record=dict(subject='528',source=source+1,dt=dt,parameters=p.tolist(),abdomenScale=g,
                objective=best[0],metrics=metrics,evaluations=calls,seconds=time.time()-start,
                inputPolicy='Recorded rect_abd_r or separately trained abdominal MAP, constant amplitude gain only',
                movementPD=False,standingPD=True,angleBoundaries=False,numberOfPhysicalCoefficients=24,
                independentValidation=False,referenceUse='Training loss and evaluation only')
            checkpoint.write_text(json.dumps(record,indent=2));(HERE/f'strict_{label}_{dt}.json').write_text(json.dumps(record,indent=2))
            print(label,'STAGE',dt,'score',round(best[0],3),'gain',round(g,3),metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
    return record

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--iterations',type=int,default=50);parser.add_argument('--sources',nargs='+',type=int,default=[1,2]);args=parser.parse_args()
    with ProcessPoolExecutor(max_workers=len(args.sources)) as pool:
        for f in as_completed([pool.submit(task,s-1,args.iterations) for s in args.sources]):f.result()
