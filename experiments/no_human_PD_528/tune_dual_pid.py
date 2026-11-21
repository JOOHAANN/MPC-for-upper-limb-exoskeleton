"""Calibrate the external angle-PID branch on all nine recorded conditions."""
from pathlib import Path
import sys,ctypes,subprocess,json
import numpy as np
from scipy.io import loadmat
from scipy.optimize import differential_evolution
from concurrent.futures import ThreadPoolExecutor
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent.parent
sys.path.insert(0,str(ROOT/'experiments/pid_damping_528'))
import tune_damping as old

def main():
    source=HERE/'dual_pid.c'; dll=HERE/'dual_pid.dll'
    subprocess.run(['C:/Library/CPP/mingw64-posix/bin/gcc.exe','-O3','-std=c99','-ffp-contract=off',
        '-shared','-static-libgcc',str(source),'-o',str(dll),'-lm'],check=True)
    lib=ctypes.CDLL(str(dll)); array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    lib.dual_pid.argtypes=[array]*10+[ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_int,
        ctypes.c_double,ctypes.c_double,array,array]
    cases=old.previous.prepare()
    for c in cases:
        c['r']=c['r'].copy(); c['r'][3:5]=0
        c['cfg']=c['cfg'].copy();c['cfg'][1:4]=0
        c['baseline']=old.previous.roll(c,np.array([0.,0.,0.,.1,1.]))
        c['threshold'][:]=np.nan; limits=[]
        for seg in c['q']['segments']:
            ix=slice(int(seg['startIndex'])-1,int(seg['endIndex']))
            limit=.75*np.max(c['baseline'][ix,3]);limits.append(limit);c['threshold'][ix]=limit
        c['threshold'][np.isnan(c['threshold'])]=min(limits)
    pid=loadmat(ROOT/'Demo/NMPC/data/pid_damped_528.mat',simplify_cells=True)['pid']
    tasks=loadmat(ROOT/'Demo/NMPC/data/no_human_PD_task_528.mat',simplify_cells=True)['tasks']
    for c in cases:
        planned=tasks[c['group']-1]
        c['qref']=np.ascontiguousarray(planned['theta']);c['vref']=np.ascontiguousarray(planned['velocity'])
        c['threshold']=np.ascontiguousarray(planned['thresholds'][:,c['source']-1])
    fixed=np.array([pid[k] for k in ('Kp','Ki','Kd','derivativeTime')])
    def rollout(c,g):
        out=np.zeros((len(c['mass']),27));q=c['q']
        failure=lib.dual_pid(c['p'],c['r'],c['cfg'],c['a'],c['mass'],c['waiting'],c['threshold'],
            c['initial'],c['qref'],c['vref'],len(out),q['dt'],10,2,80.,200.,g,out)
        if failure: raise RuntimeError(failure)
        return out
    # The C search must reproduce the MATLAB smoke run before tuning.
    smoke=loadmat(ROOT/'Demo/NMPC/output/no_human_PD_528/group1_source1.mat',simplify_cells=True)
    c=cases[0].copy();n=len(smoke['inputs'])
    for key in ('mass','waiting','a','threshold','qref','vref'):
        c[key]=np.ascontiguousarray(c[key][:n])
    c['threshold']=np.ascontiguousarray(smoke['pair']['threshold'])
    c['qref']=np.ascontiguousarray(smoke['task']['theta']);c['vref']=np.ascontiguousarray(smoke['task']['velocity'])
    c['a']=np.ascontiguousarray(smoke['inputs'])
    parity=np.max(np.abs(rollout(c,np.r_[fixed,100.,8.,25.])-smoke['pair']['pid']['actual']))
    assert parity<.001,parity
    print('PARITY',parity,flush=True)
    pool=ThreadPoolExecutor(max_workers=6);best=[np.inf,None];evaluations=[]
    def objective(x):
        g=np.r_[fixed,10**x]
        try: outputs=list(pool.map(lambda c:rollout(c,g),cases))
        except RuntimeError:return 1e9
        rows=[]
        for c,out in zip(cases,outputs):
            for seg in c['q']['segments']:
                ix=slice(int(seg['startIndex'])-1,int(seg['endIndex']))
                rows.append(dict(group=c['group'],source=c['source'],load=seg['loadKg'],speed=seg['speedCode'],
                    thetaRMSEDeg=float(np.rad2deg(np.sqrt(np.mean((out[ix,0]-c['qref'][ix])**2)))),
                    peakExcessNm=float(max(0,np.max(out[ix,3]-c['threshold'][ix]))),
                    excessAreaNms=float(np.sum(np.maximum(out[ix,3]-c['threshold'][ix],0))*.01),
                    variation=float(np.sum(np.abs(np.diff(out[ix,24]))))))
        score=np.mean([r['thetaRMSEDeg']+.08*r['peakExcessNm']+.06*r['excessAreaNms']+.002*r['variation'] for r in rows])
        evaluations.append(dict(score=float(score),gains=g.tolist()))
        if score<best[0]:
            best[:]=[score,g];print('BEST',len(evaluations),score,g[4:],flush=True)
            (HERE/'pid_selection.json').write_text(json.dumps(dict(gains=g.tolist(),score=float(score),rows=rows),indent=2))
        return score
    for g in ([100.,8.,25.],[200.,8.,45.],[300.,8.,60.]): objective(np.log10(g))
    differential_evolution(objective,[(1.8,2.8),(-1,1.3),(1.1,2.1)],maxiter=5,popsize=4,seed=528,polish=False)
    selected=json.loads((HERE/'pid_selection.json').read_text());selected['evaluations']=evaluations
    selected['calibrationCases']=18;selected['independentValidation']=False
    (HERE/'pid_selection.json').write_text(json.dumps(selected,indent=2))

if __name__=='__main__':main()
