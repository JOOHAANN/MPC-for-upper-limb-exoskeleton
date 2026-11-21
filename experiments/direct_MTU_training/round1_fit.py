"""Subject/source-specific calibration of the direct-activation dual-MTU model.

The C routine is the same source used by the readable MATLAB demonstration.
Reference signals enter this objective only, never the forward simulator.
"""
import argparse
import ctypes
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.io import loadmat, savemat
from scipy.optimize import least_squares
from scipy.signal import lfilter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEMO = ROOT / "Demo/BSE/Direct_Activation_No_Limits"
LIB = None
# Entries 22/23 encode tendon slack lengths, then decode to reference MT length.
LOW = np.array([60,.15,.45,.5,2,.4,.02,.005,5,.4,.25,.45,.345,0,.05,0,0,0,
                15,.20,.020,.035,.025,.025,2,.02,.005,.10])
HIGH = np.array([1600,1.1,1.45,45,15,3,.10,.15,1000,2,.95,1.45,.85,1.1,5,1.5,12,12,
                 45,.35,.060,.10,.35,.40,15,.10,.15,.35])


def library():
    global LIB
    if LIB is None:
        LIB = ctypes.CDLL(str(HERE / "direct_mtu.dll"))
        array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
        LIB.direct_mtu_rollout.argtypes = [array,array,ctypes.c_int,array,array,array,array,
                                           ctypes.c_double,ctypes.c_int,array]
        LIB.direct_mtu_rollout.restype = ctypes.c_int
    return LIB


def encode(p):
    p = np.array(p, dtype=float).copy()
    w = np.sin(p[27])
    p[22] -= p[20]/p[1]*np.sqrt(p[2]**2-w*w)
    p[23] -= p[21]/p[10]*np.sqrt(p[11]**2-w*w)
    return np.clip((p-LOW)/(HIGH-LOW), 1e-6, 1-1e-6)


def decode(x):
    p = LOW+np.asarray(x)*(HIGH-LOW)
    w = np.sin(p[27])
    p[22] += p[20]/p[1]*np.sqrt(p[2]**2-w*w)
    p[23] += p[21]/p[10]*np.sqrt(p[11]**2-w*w)
    return np.ascontiguousarray(p)


def memory(a, tau, dt):
    decay = np.exp(-dt/tau)
    return lfilter([1-decay], [1,-decay], a, zi=[decay*a[0]])[0]


def as_list(value):
    return [value] if isinstance(value,dict) else list(value)


def records(subject, source, dt):
    model = as_list(subject["models"])[source]
    result = []
    for q in as_list(subject["sequences"]):
        base_dt = float(q["dt"])
        back = q["backActivation"][:,source]
        fast,slow = memory(back,.1,base_dt),memory(back,.35,base_dt)
        if source == 0:
            front = q["abdomenActivation"]
        else:
            mass,load = q["mass"]/15,q["nominalLoad"]/15
            features = np.c_[back,back**2,fast,slow,back-slow,mass,load,
                             back*load,back*mass,(back-fast)*load]
            g = model["generator"]
            design = np.c_[np.ones(len(back)),(features-g["mu"])/g["scale"]]
            front = 1/(1+np.exp(-np.clip(design@g["beta"],-30,30)))
        stride = round(dt/base_dt)
        ix = np.arange(0,len(back),stride)
        time_grid = q["time"][ix]
        record = dict(back=back[ix],front=front[ix],rise=np.maximum(fast-slow,0)[ix],
                      fall=np.maximum(slow-fast,0)[ix],mass=np.ascontiguousarray(q["mass"][ix]),
                      waiting=np.ascontiguousarray(q["waiting"][ix],dtype=float),
                      initial=np.ascontiguousarray(q["initialState"],dtype=float),
                      cfg=np.ascontiguousarray(np.r_[q["standingAngle"],subject["balanceController"]]),
                      dt=dt,theta=q["referenceTheta"][ix],torque=q["referenceTorque"][ix],
                      torque_valid=q["torqueValid"][ix].astype(bool),cases=[])
        for seg in as_list(q["segments"]):
            lo = q["time"][int(seg["startIndex"])-1]
            hi = q["time"][int(seg["endIndex"])-1]
            available = np.flatnonzero((time_grid>=lo)&(time_grid<=hi))
            chosen = available[np.unique(np.round(np.linspace(0,len(available)-1,64)).astype(int))]
            tq = chosen[record["torque_valid"][chosen]&np.isfinite(record["torque"][chosen])]
            record["cases"].append((chosen,tq))
        result.append(record)
    assert sum(len(r["cases"]) for r in result)==9
    return result


def rollout(p, r):
    recruitment = p[14]*r["front"]+p[15]*r["back"]+p[16]*r["rise"]+p[17]*r["fall"]
    inputs = np.ascontiguousarray(np.c_[r["back"],-np.expm1(-np.maximum(recruitment,0))])
    output = np.full((len(inputs),27),np.nan)
    # Coarse search uses 2-ms internal steps. Final calibration uses 1 ms.
    inner = max(1,round(r["dt"]/(.002 if r["dt"]>.01 else .001)))
    last = library().direct_mtu_rollout(p,r["cfg"],len(inputs),inputs,r["mass"],r["waiting"],
                                      r["initial"],r["dt"],inner,output)
    if last<len(inputs): output[last:]=np.nan
    return output


def correlation(a,b):
    if len(a)<3 or np.std(a)<1e-10 or np.std(b)<1e-10: return 0.
    return float(np.corrcoef(a,b)[0,1])


def evaluate(p, rr, with_residual=True):
    residual=[]; qr=[];tr=[];qe=[];te=[];ratios=[];complete=True
    for r in rr:
        out=rollout(p,r)
        for ix,valid_torque in r["cases"]:
            pred_q=out[ix,0];ref_q=r["theta"][ix];scale=np.sqrt(len(ix))
            finite=np.all(np.isfinite(pred_q))
            if not finite: complete=False
            pred_q=np.nan_to_num(pred_q,nan=10,posinf=10,neginf=-10)
            residual.extend(2*(pred_q-ref_q)/(.3*scale))
            z=(pred_q-pred_q.mean())/max(np.std(pred_q),.12)
            zref=(ref_q-ref_q.mean())/max(np.std(ref_q),.12)
            residual.extend(1.5*(z-zref)/scale)
            qr.append(correlation(pred_q,ref_q));qe.append(float(np.sqrt(np.mean((pred_q-ref_q)**2))))
            if len(valid_torque):
                pred_t=np.nan_to_num(out[valid_torque,2],nan=3000,posinf=3000,neginf=-3000)
                ref_t=r["torque"][valid_torque];st=np.sqrt(len(pred_t))
                residual.extend((pred_t-ref_t)/(50*st))
                z=(pred_t-pred_t.mean())/max(np.std(pred_t),15)
                zref=(ref_t-ref_t.mean())/max(np.std(ref_t),15)
                residual.extend(2*(z-zref)/st)
                tr.append(correlation(pred_t,ref_t));te.append(float(np.sqrt(np.mean((pred_t-ref_t)**2))))
            back=np.nan_to_num(out[ix,3],nan=1000);front=np.nan_to_num(out[ix,4],nan=1000)
            ratio=np.mean(front)/max(np.mean(back),1)
            ratios.append(float(ratio))
            # Broad engineering preference: abdominal torque should not
            # dominate the back torque. This is not an observed muscle ratio.
            residual.append(2*max(ratio-.4,0))
            residual.append(.1*max(np.max(front)/max(np.max(back),1)-.8,0))
    # Wide capacity priors keep net-torque cancellation from hiding large forces.
    residual.extend([.08*np.log((p[0]/p[20])/6000),.08*np.log((p[8]/p[21])/1500)])
    metrics=dict(thetaR=float(np.mean(qr)),torqueR=float(np.mean(tr)),
                 thetaRMSEDeg=float(np.rad2deg(np.mean(qe))),torqueRMSE=float(np.mean(te)),
                 abdomenBackRatio=float(np.mean(ratios)),complete=bool(complete))
    return np.asarray(residual),metrics


def seeds(subject,source):
    p=as_list(subject["models"])[source]["parameters"]
    choices=[encode(p)]
    for back,front,damp in [(.3,.3,25),(.5,.2,30),(.2,.1,35),(.7,.5,20)]:
        trial=p.copy();trial[0]*=back;trial[8]*=front;trial[3]=damp
        choices.append(encode(trial))
    neutral=p.copy();neutral[[0,8,3,5,9,7,26]]=[350,60,30,1.2,1.2,.04,.04]
    neutral[[14,15,16,17]]=[1,0,0,0]
    choices.append(encode(neutral))
    old=loadmat(ROOT/'data/processed/hill_compliant/calibrated_models.mat',simplify_cells=True)
    try:
        choices.append(encode(as_list(old['models'])[source]['parameters'][:28]))
    except (KeyError,TypeError): pass
    return choices


def fit_route(job):
    name,source,quick=job
    subject=loadmat(HERE/'before_training/data'/f'subject_{name}.mat',simplify_cells=True)['subject']
    label=f'{name}_{source+1}';start=time.time();calls=0
    rr=records(subject,source,.04)
    candidates=[]
    for x in seeds(subject,source):
        residual,metrics=evaluate(decode(x),rr)
        candidates.append((float(residual@residual),x,metrics))
    candidates.sort(key=lambda v:v[0])
    print(label,'initial',candidates[0][0],candidates[0][2],flush=True)
    best=candidates[0]
    # Full 28-variable optimization; no movement feedback or boundary variables.
    for seed_index,(_,x,_) in enumerate(candidates[:1 if quick else 2]):
        stage_best=[np.inf,x.copy()]
        def objective(z):
            nonlocal calls
            calls+=1
            residual,_=evaluate(decode(z),rr)
            score=float(residual@residual)
            if score<stage_best[0]:stage_best[:]=[score,z.copy()]
            return residual
        result=least_squares(objective,x,bounds=(np.zeros(28),np.ones(28)),
                             max_nfev=18 if quick else 45,diff_step=2e-4,
                             ftol=2e-4,xtol=2e-4,gtol=1e-4,x_scale='jac')
        score,x=stage_best;_,metrics=evaluate(decode(x),rr)
        print(label,'coarse',seed_index,score,metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
        if score<best[0]:best=(score,x,metrics)
    stages=[(.02,15 if quick else 40),(.01,10 if quick else 35)]
    x=best[1]
    for dt,iterations in stages:
        rr=records(subject,source,dt);stage_best=[np.inf,x.copy()]
        result=least_squares(objective,x,bounds=(np.zeros(28),np.ones(28)),
                             max_nfev=iterations,diff_step=1e-4,
                             ftol=1e-4,xtol=1e-4,gtol=1e-4,x_scale='jac')
        score,x=stage_best;residual,metrics=evaluate(decode(x),rr)
        record=dict(subject=name,source=source+1,dt=dt,objective=score,metrics=metrics,
                    parameters=decode(x).tolist(),evaluations=calls,seconds=time.time()-start,
                    status=int(result.status),message=str(result.message))
        (HERE/f'checkpoint_{label}_{dt:.2f}.json').write_text(json.dumps(record,indent=2))
        print(label,'refined',dt,score,metrics,'calls',calls,'seconds',round(time.time()-start),flush=True)
    savemat(HERE/f'fit_{label}.mat',dict(parameters=decode(x),source=source+1,subject=name,
                                      objective=score,evaluations=calls,elapsed=time.time()-start))
    (HERE/f'fit_{label}.json').write_text(json.dumps(record,indent=2))
    return record


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--subjects',nargs='+',default=['528','617','725','730'])
    parser.add_argument('--sources',nargs='+',type=int,default=[1,2]);parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--quick',action='store_true');parser.add_argument('--benchmark',action='store_true')
    args=parser.parse_args()
    if args.benchmark:
        subject=loadmat(HERE/'before_training/data/subject_528.mat',simplify_cells=True)['subject']
        for dt in [.04,.02,.01]:
            rr=records(subject,0,dt);p=decode(seeds(subject,0)[1]);tick=time.time()
            residual,metrics=evaluate(p,rr);print(dt,time.time()-tick,float(residual@residual),metrics,flush=True)
        return
    jobs=[(n,s-1,args.quick) for n in args.subjects for s in args.sources]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures=[executor.submit(fit_route,job) for job in jobs]
        for future in as_completed(futures):
            record=future.result();print('DONE',record['subject'],record['source'],record['metrics'],flush=True)


if __name__=='__main__':main()
