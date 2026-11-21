"""Tune a single reactive PID on all 18 current-model source/condition cases."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import ctypes,json,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.io import loadmat
from scipy.optimize import minimize
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
data=loadmat(ROOT/'Demo/BSE/Direct_Activation_No_Limits/data/subject_528.mat',simplify_cells=True)['subject']
baseline=loadmat(ROOT/'Demo/BSE/Direct_Activation_No_Limits/output/BSE_results.mat',simplify_cells=True)['results']['subjects']
if not isinstance(baseline,dict):baseline=baseline[0]
arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
lib=ctypes.CDLL(str(HERE/'pid_tuning_filtered.dll'))
lib.current_pid_run.argtypes=[arr,arr,arr,arr,arr,ctypes.c_int,arr,arr,arr,arr,arr,ctypes.c_double,ctypes.c_int,arr]
lib.current_pid_run.restype=ctypes.c_int
pd=np.array([5.,3.,.005,.3,-.001/np.log(.92)]);actuator=np.array([80.,200.]);records=[]
for group,q in enumerate(data['sequences']):
    for source in range(2):
        p=np.ascontiguousarray(data['models'][source]['parameters']);back=q['backActivation'][:,source]
        front=q['abdomenActivation' if source==0 else 'abdomenMAPActivation'];inputs=np.ascontiguousarray(np.clip(np.c_[back,front],0,1))
        base=baseline['groups'][group]['runs'][source]['output'];threshold=np.zeros(len(back));cases=[];limits=[]
        for seg in q['segments']:
            ix=np.arange(int(seg['startIndex'])-1,int(seg['endIndex']));limit=.75*base[ix,3].max();limits.append(limit);threshold[ix]=limit;cases.append(ix)
        threshold[q['waiting']>.5]=min(limits)
        records.append(dict(p=p,cfg=np.ascontiguousarray(np.r_[q['standingAngle'],data['balanceController']]),
            input=inputs,mass=np.ascontiguousarray(q['mass']),waiting=np.ascontiguousarray(q['waiting'],dtype=float),
            initial=np.ascontiguousarray(q['initialState']),threshold=threshold,baseline=base,cases=cases,dt=float(q['dt'])))

baselineArea=sum(float(np.maximum(r['baseline'][ix,3]-r['threshold'][ix],0).sum()*r['dt']) for r in records for ix in r['cases'])

def evaluate(pid,details=False):
    pid=np.ascontiguousarray(pid);scores=[];qerr=[];excess=[];area=[]
    for r in records:
        out=np.empty((len(r['input']),27));n=len(out)
        done=lib.current_pid_run(r['p'],r['cfg'],pd,actuator,pid,n,r['input'],r['mass'],r['waiting'],r['initial'],r['threshold'],r['dt'],10,out)
        if done!=n:return (1e7,{}) if details else 1e7
        for ix in r['cases']:
            target=np.minimum(r['baseline'][ix,3],r['threshold'][ix]);t=out[ix,3]-target
            q=out[ix,0]-r['baseline'][ix,0];u=out[ix,24];e=np.maximum(out[ix,3]-r['threshold'][ix],0)
            cost=np.mean((t/30)**2)+2*np.mean((q/.1)**2)+(e.max()/30)**2+.05*np.sum(np.abs(np.diff(u)))/80
            scores.append(float(cost));qerr.append(float(np.rad2deg(np.sqrt(np.mean(q*q)))));excess.append(float(e.max()));area.append(float(e.sum()*r['dt']))
    value=float(np.mean(scores))+200*max(sum(area)/baselineArea-.6,0)**2
    info=dict(meanThetaRMSEDeg=float(np.mean(qerr)),maximumThetaRMSEDeg=max(qerr),maximumExcessNm=max(excess),totalExcessAreaNms=sum(area),excessAreaRatio=sum(area)/baselineArea)
    return (value,info) if details else value

def decode(z):return 10**np.asarray(z)
def main():
    start=time.time();rng=np.random.default_rng(52892);low=np.log10([.05,.05,.0001,.03,.03]);high=np.log10([3,40,.2,.8,.4])
    candidates=[np.array([.4,3,.008,.2,.15]),np.array([.4132547,6.35203,.0120483,.467284,.15])]
    candidates.extend([np.array([kp,ki,.01,.3,release]) for kp in [.2,.5,1.] for ki in [.5,3.,10.] for release in [.08,.18]])
    candidates.extend(decode(rng.uniform(low,high,(48,5))))
    with ThreadPoolExecutor(max_workers=8) as pool:values=list(pool.map(evaluate,candidates))
    order=np.argsort(values);best=[float(values[order[0]]),candidates[order[0]]]
    for i in order[:5]:print('CANDIDATE',values[i],candidates[i].tolist(),evaluate(candidates[i],True)[1],flush=True)
    def objective(z):
        pid=decode(z);value=evaluate(pid)
        if value<best[0]:best[:]=[value,pid.copy()]
        return value
    initial=np.maximum(best[1],10**low);
    result=minimize(objective,np.log10(initial),method='Nelder-Mead',bounds=list(zip(low,high)),options={'maxiter':40,'xatol':.015,'fatol':.01})
    value,info=evaluate(best[1],True);record=dict(pid=dict(zip(['Kp','Ki','Kd','derivativeTime','releaseTime'],best[1].tolist())),objective=value,metrics=info,
        objectiveDefinition='mean((torque_error/30)^2) + 2*mean((theta_deviation/0.1)^2) + (peak_excess/30)^2 + 0.05*assistance_total_variation/80',
        baselineExcessAreaNms=baselineArea,minimumDesiredExcessAreaReductionFraction=.4,
        releasePolicy='Filter command withdrawal only; onset is unfiltered. Set command to zero below 0.05 N m. Hardware withdrawal remains instantaneous.',
        conditions=18,oneTuningAcrossBothSources=True,matchedActuator=[80,200],seconds=time.time()-start)
    (HERE/'pid_fit.json').write_text(json.dumps(record,indent=2));print('FINAL',json.dumps(record),flush=True)
if __name__=='__main__':main()
