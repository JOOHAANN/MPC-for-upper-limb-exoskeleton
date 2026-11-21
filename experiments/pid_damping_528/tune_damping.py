"""Retune the causal PID for damping under the shared asymmetric actuator."""
from pathlib import Path
import importlib.util, json, time, argparse
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.io import loadmat
from scipy.signal import butter, sosfiltfilt, find_peaks
from scipy.optimize import differential_evolution

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
spec=importlib.util.spec_from_file_location('previous',HERE.parent/'controller_comparison_528/tune_pid.py')
previous=importlib.util.module_from_spec(spec); spec.loader.exec_module(previous)
previous.HERE=HERE
previous.ROOT=ROOT
previous.DEMO=ROOT/'Demo'
OLD=np.array([.402297027908,60.1251089985,.0389523012315,.0659541494825,1.])
FILTER=butter(3,3,fs=100,btype='highpass',output='sos')

def measures(c,out):
    rows=previous.metrics(c,out)
    for row,seg in zip(rows,c['q']['segments']):
        ix=slice(int(seg['startIndex'])-1,int(seg['endIndex']))
        u=out[ix,24]; back=out[ix,3]
        # High-frequency energy is a ripple diagnostic, not a stability proof.
        row['ripple']=float(np.sqrt(np.mean(sosfiltfilt(FILTER,u)**2)))
        row['back_ripple']=float(np.sqrt(np.mean(sosfiltfilt(FILTER,back)**2)))
        row['variation']=float(np.sum(np.abs(np.diff(u))))
        row['reversals']=int(len(find_peaks(u,prominence=2)[0])+len(find_peaks(-u,prominence=2)[0]))
    return rows

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--iterations',type=int,default=14)
    args=parser.parse_args(); start=time.monotonic(); cases=previous.prepare()
    saved=loadmat(ROOT/'Demo/NMPC/output/pid_nmpc_528/threshold_tracking/pid_nmpc_comparison_528.mat',simplify_cells=True)['results']
    global OLD
    OLD=np.array([saved['pid'][k] for k in ('Kp','Ki','Kd','derivativeTime')]+[1.])
    oldrows=[]
    for c in cases:
        pair=saved['groups'][c['group']-1]['runs'][c['source']-1]
        c['threshold']=np.ascontiguousarray(pair['threshold'])
        c['a']=np.ascontiguousarray(pair['nmpc']['inputs'])
        c['baseline']=pair['baseline']
        actual=pair['pid']['actual']
        oldrows+=measures(c,actual)
    (HERE/'before_metrics.json').write_text(json.dumps(oldrows,indent=2))
    print('OLD',json.dumps({key:float(np.mean([r[key] for r in oldrows])) for key in ('peak','area','ripple','back_ripple','variation','reversals')}),flush=True)
    best=[float('inf'),None]; evaluations=[]
    pool=ThreadPoolExecutor(max_workers=6)
    def objective(x):
        g=np.r_[10**x,1.]
        try: rows=[r for block in pool.map(lambda c: measures(c,previous.roll(c,g)),cases) for r in block]
        except RuntimeError: return 1e8
        # Preserve useful assistance while prioritizing smooth, damped action.
        score=np.mean([.12*r['peak']+.2*r['area']+1.8*r['ripple']+.012*r['variation']+
                       .15*r['reversals']+8*r['theta'] for r in rows])
        score+=.8*max(r['ripple'] for r in rows)
        evaluations.append(dict(score=float(score),gains=g.tolist()))
        if score<best[0]:
            best[:]=[score,g.copy()]
            print('BEST',len(evaluations),round(time.monotonic()-start,1),round(score,5),g,flush=True)
            (HERE/'selection.json').write_text(json.dumps(dict(gains=g.tolist(),score=float(score),rows=rows),indent=2))
        return score
    for gains in ([.5,8,.02,.2],[.5,5,.01,.25],[.4,8,.005,.25],[.4,2,.02,.1]):
        objective(np.log10(gains))
    result=differential_evolution(objective,[(-.6,.2),(0,1.3),(-3,-1.3),(-1,-.3)],
        maxiter=args.iterations,popsize=5,seed=528,polish=False,tol=.003)
    selected=json.loads((HERE/'selection.json').read_text()); selected['evaluations']=evaluations
    selected['oldGains']=OLD.tolist(); selected['tuningConditions']=18
    selected['objective']='0.12 peak + 0.2 area + 1.8 HF ripple + .012 TV + .15 reversals + 8 theta RMSE; plus .8 worst ripple'
    (HERE/'selection.json').write_text(json.dumps(selected,indent=2))
    print('FINAL',selected['gains'],selected['score'],flush=True)

if __name__=='__main__': main()
