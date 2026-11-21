"""Train a genuine abdominal MAP: NMF4 -> Gaussian3 -> condition weights.

Only experiment-derived rect_abd_r activation and activation-based phase
windows are used. No IK, reference torque or BSE output is read for fitting.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.io import loadmat,savemat
from scipy.ndimage import gaussian_filter1d,binary_closing
from scipy.optimize import least_squares,nnls
from scipy.signal import lfilter,find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
DEMO=ROOT/'Demo/BSE/Direct_Activation_No_Limits'
RAW=ROOT.parent/'YOUHAN LI code and project description/Data/240528Pilot'
OUT=HERE/'abdominal_MAP_528';OUT.mkdir(exist_ok=True)
PHASE=np.linspace(0,1,200)
REP=['IL_L1_r','IL_L1_l','LTpT_T7_r','LTpT_T7_l','LTpL_L1_r','LTpL_L1_l']


def cells(value):return [value] if isinstance(value,dict) else list(value)


def read(path,names=None):
    with path.open() as stream:
        header=[]
        for line in stream:
            header.append(line)
            if line.strip()=='endheader':break
        labels=stream.readline().split()
    columns=[0]+[labels.index(n) for n in names] if names else [0]
    a=np.loadtxt(path,skiprows=len(header)+1,usecols=columns)
    return labels,a


def known_windows():
    subject=loadmat(HERE/'before_training/data/subject_528.mat',simplify_cells=True)['subject']
    result={}
    for q in cells(subject['sequences']):
        offset=0
        for stage,seg in enumerate(cells(q['segments'])):
            load,speed=int(seg['loadKg']),int(seg['speedCode']);base=RAW/f'dynamics_stoop_{load}kg_{speed}'
            _,time=read(base/'lmt.sto')
            first,last=(49 if stage==2 else 0),(349 if stage==1 else 299)
            ix=np.arange(int(seg['startIndex'])-1,int(seg['endIndex']))
            epoch=time[first]+q['time'][ix]-offset*.001
            result[(float(load),speed)]=(epoch[0],epoch[-1],q['abdomenActivation'][ix].copy())
            offset+=int(np.floor((time[last]-time[first])/.001+.0003))+1
    return result


def collect():
    known=known_windows();records=[]
    for token,load in [('5',5.),('75',7.5),('10',10.),('125',12.5),('15',15.)]:
        for speed in [25,35,45,55,65]:
            path=RAW/f'dynamics_stoop_{token}kg_{speed}'/'Activations.sto'
            _,a=read(path,['rect_abd_r',*REP]);a=a[np.unique(a[:,0],return_index=True)[1]]
            t=np.arange(a[0,0],a[-1,0]-.005,.01)
            front=np.interp(t,a[:,0],a[:,1]);back=np.mean(np.column_stack([np.interp(t,a[:,0],a[:,j]) for j in range(2,8)]),axis=1)
            d=np.exp(-.01/.12);front=lfilter([1-d],[1,-d],front,zi=[d*front[0]])[0]
            if (load,speed) in known:
                lo,hi,measured=known[(load,speed)];duration=hi-lo
                profile=np.interp(PHASE,np.linspace(0,1,len(measured)),measured)
                window_source='Existing experimental-activation bout; no IK/torque used'
            else:
                env=gaussian_filter1d(back,10);baseline,peak=np.percentile(env,[15,99])
                active=env>baseline+.15*(peak-baseline)
                active=binary_closing(active,structure=np.ones(71))
                ix=np.flatnonzero(active)
                if len(ix)<50 or peak-baseline<.01:raise ValueError(f'No usable activation bout: {path}')
                first=max(0,ix[0]-20);last=min(len(t)-1,ix[-1]+20)
                if t[last]-t[first]>18:
                    peaks,properties=find_peaks(env,prominence=.15*(peak-baseline),distance=60)
                    pairs=[j for j in range(len(peaks)-1) if .5<t[peaks[j+1]]-t[peaks[j]]<10]
                    if not pairs:raise ValueError(f'Cannot identify one activation pair: {path}')
                    j=max(pairs,key=lambda j:min(properties['prominences'][j:j+2]))
                    first=max(0,peaks[j]-50);last=min(len(t)-1,peaks[j+1]+50)
                lo,hi=t[first],t[last];duration=hi-lo
                profile=np.interp(PHASE,(t[first:last+1]-lo)/duration,front[first:last+1])
                window_source='Back-activation envelope only; no IK/torque used'
            assert duration>1 and np.all(np.isfinite(profile)) and profile.min()>=0
            records.append(dict(loadKg=load,speedCode=speed,cadenceBPM=60/duration,duration=duration,
                profile=profile,sourcePath=str(path),sourceSHA256=hashlib.sha256(path.read_bytes()).hexdigest(),
                sourceStartEpoch=float(lo),sourceEndEpoch=float(hi),windowSource=window_source))
    assert len(records)==25 and len(set(r['loadKg'] for r in records))==5
    return records


def gaussian(z):
    return sum(a*np.exp(-((PHASE-c)/w)**2) for a,c,w in np.asarray(z).reshape(3,3))


def features(load,cadence):return np.array([load/15,(load/15)**2,cadence/10,(cadence/10)**2])


def train(records):
    # A single global channel scale preserves between-condition amplitudes.
    amplitude=max(float(np.percentile(np.concatenate([r['profile'] for r in records]),99)),.01)
    y=np.column_stack([r['profile']/amplitude for r in records]);rng=np.random.default_rng(528)
    w=rng.uniform(.1,.8,(200,4));h=rng.uniform(.01,.3,(4,len(records)))
    for _ in range(1200):
        h*=(w.T@y)/np.maximum(w.T@w@h+.0001*h,1e-12)
        w*=(y@h.T)/np.maximum(w@h@h.T+.0001*w,1e-12)
        scale=np.maximum(np.max(w,axis=0),1e-10);w/=scale;h*=scale[:,None]
    parameters=[]
    for j in range(4):
        best=None
        for centers in [(.2,.5,.8),(.3,.45,.7),(.15,.65,.85)]:
            z=np.array([[.6,c,.15] for c in centers]).ravel()
            result=least_squares(lambda z:gaussian(z)-w[:,j],z,
                bounds=(np.tile([-2,0,.025],3),np.tile([3,1,.6],3)),max_nfev=500)
            if best is None or result.cost<best.cost:best=result
        parameters.append(best.x)
    primitives=np.maximum(0,np.column_stack([gaussian(z) for z in parameters]))
    target=np.array([nnls(primitives,r['profile']/amplitude)[0] for r in records])
    x=np.array([features(r['loadKg'],r['cadenceBPM']) for r in records]);mu=x.mean(0);scale=np.maximum(x.std(0),.1)
    design=np.c_[np.ones(len(x)),(x-mu)/scale];penalty=np.eye(5)*.15;penalty[0,0]=1e-8
    beta=np.linalg.solve(design.T@design+penalty,design.T@target)
    return dict(primitives=primitives,gaussianParameters=np.asarray(parameters),mu=mu,scale=scale,beta=beta,
        covariateMinimum=x.min(0),covariateMaximum=x.max(0),phase=PHASE,amplitudeScale=amplitude,
        targetMuscle='rect_abd_r',architecture='4 NMF temporal primitives, 3 Gaussian lobes per primitive, quadratic load/cadence weight regression',
        usesIK=False,usesTorque=False,usesBSEOutput=False,usesBackToAbdomenRecruitment=False)


def predict(model,r):
    x=np.clip(features(r['loadKg'],r['cadenceBPM']),model['covariateMinimum'],model['covariateMaximum'])
    weights=np.maximum(0,np.r_[1,(x-model['mu'])/model['scale']]@model['beta'])
    return np.clip(model['amplitudeScale']*(model['primitives']@weights),0,1)


def main():
    records=collect();hold=[r for r in records if r['loadKg'] in [7.5,12.5]]
    initial=train([r for r in records if r['loadKg'] not in [7.5,12.5]])
    holdout=[]
    for r in hold:
        pred=predict(initial,r);corr=float(np.corrcoef(pred,r['profile'])[0,1])
        holdout.append(dict(loadKg=r['loadKg'],speedCode=r['speedCode'],rmse=float(np.sqrt(np.mean((pred-r['profile'])**2))),correlation=corr))
    model=train(records);rows=[]
    fig,axes=plt.subplots(5,5,figsize=(15,12),sharex=True)
    for r,ax in zip(records,axes.ravel()):
        pred=predict(model,r);corr=float(np.corrcoef(pred,r['profile'])[0,1])
        rows.append(dict(loadKg=r['loadKg'],speedCode=r['speedCode'],cadenceBPM=r['cadenceBPM'],
            rmse=float(np.sqrt(np.mean((pred-r['profile'])**2))),correlation=corr))
        ax.plot(PHASE,r['profile'],'k--',lw=1,label='Experimental-derived abdomen')
        ax.plot(PHASE,pred,color='#7b3294',lw=1.2,label='Trained abdominal MAP')
        ax.set_title(f"{r['loadKg']:g} kg | code {r['speedCode']}",fontsize=9);ax.grid(alpha=.2)
    axes[-1,2].set_xlabel('Normalized activation-bout phase');axes[2,0].set_ylabel('Abdominal activation input')
    fig.suptitle('Subject 528 | Abdominal MAP trained from experimental rect_abd_r | All five load knots',fontsize=13)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2)
    fig.tight_layout(rect=[0,.03,1,.97]);fig.savefig(OUT/'abdominal_MAP_25_conditions.png',dpi=180);plt.close(fig)
    record_cells=np.empty((1,len(records)),dtype=object)
    for i,r in enumerate(records):record_cells[0,i]=r
    savemat(DEMO/'data/abdominal_MAP_528.mat',dict(model=model,records=record_cells),do_compression=True,long_field_names=True)
    savemat(OUT/'model.mat',dict(model=model),do_compression=True,long_field_names=True)
    with (OUT/'condition_metrics.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    audit=dict(subject='528',trainingConditions=25,loadKnots=[5,7.5,10,12.5,15],
        finalMeanRMSE=float(np.mean([r['rmse'] for r in rows])),finalMeanCorrelation=float(np.mean([r['correlation'] for r in rows])),
        heldOutLoadsBeforeFinalRefit=[7.5,12.5],heldOutMeanRMSE=float(np.mean([r['rmse'] for r in holdout])),
        heldOutMeanCorrelation=float(np.mean([r['correlation'] for r in holdout])),holdout=holdout,
        architecture=model['architecture'],usesIK=False,usesTorque=False,usesBSEOutput=False,
        provenance=[{k:v for k,v in r.items() if k!='profile'} for r in records])
    (OUT/'training_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps({k:v for k,v in audit.items() if k not in ['provenance','holdout']},indent=2))


if __name__=='__main__':main()
