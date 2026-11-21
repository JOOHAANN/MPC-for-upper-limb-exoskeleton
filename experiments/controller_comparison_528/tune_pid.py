"""Develop a fixed threshold-triggered PID without changing the frozen demo."""
from pathlib import Path
import argparse, ctypes, hashlib, json, subprocess, sys, time
import numpy as np
from scipy.io import loadmat
from scipy.signal import lfilter
from scipy.optimize import differential_evolution

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DEMO = ROOT / 'Demo'

def memory(x,tau,dt):
    d=np.exp(-dt/tau)
    return lfilter([1-d],[1,-d],x,zi=[d*x[0]])[0]

def inputs(q,m,source):
    back=q['backActivation'][:,source]
    if source==0:
        front=q['abdomenActivation']
    else:
        fast=memory(back,.10,q['dt']); slow=memory(back,.35,q['dt'])
        load=q['nominalLoad']/15; mass=q['mass']/15
        f=np.column_stack([back,back**2,fast,slow,back-slow,mass,load,
                          back*load,back*mass,(back-fast)*load])
        g=m['generator']; design=np.column_stack([np.ones(len(back)),(f-g['mu'])/g['scale']])
        logit=np.clip(design@g['beta'],-30,30); front=1/(1+np.exp(-logit))
    p=m['parameters']; change=memory(back,.1,q['dt'])-memory(back,.35,q['dt'])
    recruitment=p[14]*front+p[15]*back+p[16]*np.maximum(change,0)+p[17]*np.maximum(-change,0)
    return np.ascontiguousarray(np.column_stack([back,-np.expm1(-np.maximum(recruitment,0))]))

def prepare():
    source=HERE/'pid_tuning_core.c'; dll=HERE/'pid_tuning_core.dll'
    if not dll.exists() or dll.stat().st_mtime<source.stat().st_mtime:
        subprocess.run(['C:/Library/CPP/mingw64-posix/bin/gcc.exe','-O3','-std=c99',
                        '-ffp-contract=off','-shared','-static-libgcc',str(source),'-o',str(dll),'-lm'],check=True)
    lib=ctypes.CDLL(str(dll)); array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    lib.compare_pid_rollout.argtypes=([array]*8+[ctypes.c_int,ctypes.c_double,ctypes.c_int,
                                          ctypes.c_int,ctypes.c_double,ctypes.c_double,array,array])
    subject=loadmat(DEMO/'BSE/data/subject_528.mat',simplify_cells=True)['subject']
    protected={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
               for base in (DEMO/'BSE',DEMO/'NMPC') for p in base.rglob('*')
               if p.is_file() and 'output' not in p.parts and p.suffix in ('.m','.mat','.csv','.c','.h')}
    snapshot=HERE/'original_file_hashes.json'
    if not snapshot.exists(): snapshot.write_text(json.dumps(protected,indent=2))
    cases=[]
    for group,q in enumerate(subject['sequences']):
        for src,m in enumerate(subject['models']):
            p=np.ascontiguousarray(m['parameters']); r=np.ascontiguousarray(m['neuralParameters'])
            cfg=np.ascontiguousarray(np.r_[q['standingAngle'],subject['balanceController']])
            mass=np.ascontiguousarray(q['mass'],dtype=float); waiting=np.ascontiguousarray(q['waiting'],dtype=float)
            initial=np.ascontiguousarray(q['initialState']); a=inputs(q,m,src)
            c=dict(group=group+1,source=src+1,q=q,p=p,r=r,cfg=cfg,mass=mass,
                   waiting=waiting,initial=initial,a=a,lib=lib)
            c['threshold']=np.ones(len(mass))*100
            baseline=roll(c,np.array([0.,0.,0.,.04,1.]))
            limits=[]
            for seg in q['segments']:
                idx=slice(int(seg['startIndex'])-1,int(seg['endIndex']))
                limit=.75*np.max(baseline[idx,3]); limits.append(limit); c['threshold'][idx]=limit
            c['threshold'][waiting>.5]=min(limits); c['baseline']=baseline
            cases.append(c)
    return cases

def roll(c,gains):
    out=np.zeros((len(c['mass']),27)); q=c['q']
    flag=c['lib'].compare_pid_rollout(c['p'],c['r'],c['cfg'],c['a'],c['mass'],c['waiting'],
                  c['threshold'],c['initial'],len(out),q['dt'],10,2,80.,200.,np.ascontiguousarray(gains),out)
    if flag: raise RuntimeError(f'Hill failure at {flag}')
    return out

def metrics(c,out):
    rows=[]; q=c['q']; dt=q['dt']
    for seg in q['segments']:
        idx=np.arange(int(seg['startIndex'])-1,int(seg['endIndex']))
        # Include the complete release tail before the next action.
        start=max(0,idx[0]-100); end=min(len(out),idx[-1]+151); full=slice(start,end)
        threshold=c['threshold'][idx[0]]
        excess=np.maximum(out[idx,3]-threshold,0)
        rows.append(dict(group=c['group'],source=c['source'],load=seg['loadKg'],speed=seg['speedCode'],
             limit=float(threshold),peak=float(np.max(excess)),area=float(np.sum(excess)*dt),
             theta=float(np.sqrt(np.mean((out[idx,0]-c['baseline'][idx,0])**2))),
             effort=float(np.sum(out[full,24]**2)*dt),u_peak=float(np.max(out[full,24])),
             mean_under=float(np.mean(np.maximum(np.minimum(c['baseline'][idx,3],threshold)-out[idx,3],0)))))
    return rows

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--optimize',action='store_true')
    parser.add_argument('--iterations',type=int,default=15); args=parser.parse_args()
    cases=prepare(); start=time.monotonic()
    gains=np.array([2.,20.,.05,.04,1.])
    print('INITIAL',json.dumps([r for c in cases for r in metrics(c,roll(c,gains))]),flush=True)
    print('BASELINE SIX-ROLLOUT SECONDS',time.monotonic()-start,flush=True)
    if not args.optimize: return
    # Tune only group 1 (three loads and both activation sources); groups 2/3 are held out.
    train=[c for c in cases if c['group']==1]
    evaluations=[]; best=[np.inf,None]
    def objective(x):
        g=np.r_[10**x[:3],10**x[3],1.]
        try: rows=[r for c in train for r in metrics(c,roll(c,g))]
        except RuntimeError: return 1e8
        score=np.mean([r['area']/10+.30*r['peak']+8*r['theta']+.0001*r['effort'] for r in rows])
        evaluations.append(dict(score=float(score),gains=g.tolist()))
        if score<best[0]:
            best[:]=[score,g.copy()]
            print('BEST',len(evaluations),round(time.monotonic()-start,1),score,g,flush=True)
            (HERE/'pid_selection.json').write_text(json.dumps(dict(gains=g.tolist(),score=score,
                tuningGroups=[1],heldOutGroups=[2,3],objective='mean(excess area/10 + .30*peak + 8*theta RMSE + .0001*effort)'),indent=2))
        return score
    result=differential_evolution(objective,[(-.6,1.1),(-.3,2.),(-3.,-.3),(-2.,-.5)],
                                  maxiter=args.iterations,popsize=5,seed=528,polish=False,tol=.002)
    gains=np.r_[10**result.x[:3],10**result.x[3],1.]
    all_rows=[r for c in cases for r in metrics(c,roll(c,gains))]
    (HERE/'pid_selection.json').write_text(json.dumps(dict(gains=gains.tolist(),score=float(result.fun),
        tuningGroups=[1],heldOutGroups=[2,3],rows=all_rows,evaluations=evaluations),indent=2))
    print('FINAL',json.dumps(all_rows),flush=True)

if __name__=='__main__': main()
