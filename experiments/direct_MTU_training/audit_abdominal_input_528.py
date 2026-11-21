"""Read-only reconstruction of the recorded abdominal input and its mixing."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.io import loadmat,savemat
from scipy.signal import lfilter
import fit_direct_MTU as fit

OUT=fit.HERE/'mixed_abdominal_input_audit';OUT.mkdir(parents=True,exist_ok=True)
RAW=fit.ROOT.parent/'YOUHAN LI code and project description/Data/240528Pilot'
CANDIDATE=fit.HERE/'mixed_abdomen_candidate'
subject=loadmat(CANDIDATE/'data/subject_528.mat',simplify_cells=True)['subject']
replay=loadmat(CANDIDATE/'output/BSE_results.mat',simplify_cells=True)['results']
result=next(s for s in fit.as_list(replay['subjects']) if '528' in s['id'])
rows=[];plots=[];totals=[np.zeros(4),np.zeros(4)]


def table(path,columns):
    with path.open() as stream:
        lines=[]
        for line in stream:
            lines.append(line)
            if line.strip()=='endheader':break
        names=stream.readline().split()
    values=np.loadtxt(path,skiprows=len(lines)+1,usecols=columns)
    return names,values


for group,q in enumerate(fit.as_list(subject['sequences'])):
    count=len(q['time']);samples=np.round(q['time']/.001).astype(int)
    recorded=np.zeros(count);offset=0;raw_ranges=[]
    for stage,seg in enumerate(fit.as_list(q['segments'])):
        load=int(seg['loadKg']);speed=int(seg['speedCode'])
        base=RAW/f'dynamics_stoop_{load}kg_{speed}'
        labels,activation=table(base/'Activations.sto',(0,1))
        _,lmt=table(base/'lmt.sto',(0,))
        assert labels[1]=='rect_abd_r'
        first,last=(49 if stage==2 else 0),(349 if stage==1 else 299)
        n=int(np.floor((lmt[last]-lmt[first])/.001+.0003))+1
        ix=np.flatnonzero((samples>=offset)&(samples<offset+n))
        query=lmt[first]+(samples[ix]-offset)*.001
        recorded[ix]=np.interp(query,activation[:,0],activation[:,1])
        active=np.arange(int(seg['startIndex'])-1,int(seg['endIndex']))
        raw_ranges.append(dict(loadKg=load,speedCode=speed,source=str(base/'Activations.sto'),
            sourceSHA256=hashlib.sha256((base/'Activations.sto').read_bytes()).hexdigest(),
            channel='rect_abd_r',sourceSamples=len(activation),
            interpolationOutsideSourceRows=int(np.sum((query<activation[0,0])|(query>activation[-1,0])))))
        offset+=n
    decay=np.exp(-q['dt']/.12)
    reconstructed=lfilter([1-decay],[1,-decay],recorded,zi=[decay*recorded[0]])[0]
    error=float(np.max(np.abs(reconstructed-q['abdomenActivation'])))
    assert error<1e-12,(group,error)
    for entry in raw_ranges:entry['cachedReconstructionMaximumError']=error;rows.append(entry)
    for source,model in enumerate(fit.as_list(subject['models'])):
        p=model['parameters'];back=q['backActivation'][:,source]
        fast=fit.memory(back,.1,q['dt']);slow=fit.memory(back,.35,q['dt']);change=fast-slow
        if source==0:base_front=q['abdomenActivation']
        else:
            mass=q['mass']/15;load=q['nominalLoad']/15
            features=np.c_[back,back**2,fast,slow,back-slow,mass,load,back*load,back*mass,(back-fast)*load]
            g=model['generator'];logit=np.c_[np.ones(len(back)),(features-g['mu'])/g['scale']]@g['beta']
            base_front=1/(1+np.exp(-np.clip(logit,-30,30)))
        terms=np.c_[p[14]*base_front,p[15]*back,p[16]*np.maximum(change,0),p[17]*np.maximum(-change,0)]
        prepared=-np.expm1(-np.maximum(terms.sum(axis=1),0));moving=q['waiting']<.5
        totals[source]+=terms[moving].sum(axis=0)
        run=fit.as_list(fit.as_list(result['groups'])[group]['runs'])[source]
        assert np.max(np.abs(prepared-run['inputs'][:,1]))<1e-12
        for stage,seg in enumerate(fit.as_list(q['segments'])):
            ix=np.arange(int(seg['startIndex'])-1,int(seg['endIndex']))
            if group==0:
                plots.append(dict(source=source+1,loadKg=seg['loadKg'],speedCode=seg['speedCode'],
                    time=q['time'][ix]-q['time'][ix[0]],recordedInterpolated=recorded[ix],
                    experimentalStored=q['abdomenActivation'][ix],sourceBase=base_front[ix],
                    preparedCommand=prepared[ix],actualActivation=run['output'][ix,14]))

shares={}
for source,values in enumerate(totals):
    shares['Experimental' if source==0 else 'MAP']=dict(zip(
        ['abdominalBase','backSignal','backRisingHistory','backFallingHistory'],(values/values.sum()).tolist()))
report=dict(subject='528',reconstructionPassed=True,experimentalChannel='rect_abd_r',
    auditedVersion='Archived mixed-abdomen candidate, not the strict-input release',
    modelPackageSHA256=hashlib.sha256((CANDIDATE/'data/subject_528.mat').read_bytes()).hexdigest(),
    replaySHA256=hashlib.sha256((CANDIDATE/'output/BSE_results.mat').read_bytes()).hexdigest(),
    recruitmentCoefficients=[m['parameters'][14:18].tolist() for m in fit.as_list(subject['models'])],
    sourceFile='Activations.sto in each original experimental condition folder',
    preprocessing='Timestamp interpolation to the archived 100-Hz sequence, then a causal 0.12-s low-pass filter.',
    actualInputPolicy='Both branches additionally mix abdominal base, back signal and back-history terms before nonlinear saturation.',
    mapAbdomenPolicy='Auxiliary logistic regression fitted to experimental abdominal labels, not a directly exported abdominal MAP channel.',
    precursorContributionFractions=shares,
    interpretation='Fractions refer to the nonnegative linear recruitment precursor, not actual activation or muscle-torque contribution.',
    records=rows)
(OUT/'abdominal_provenance_528.json').write_text(json.dumps(report,indent=2))
savemat(OUT/'abdominal_audit_curves.mat',dict(curves=np.array(plots,dtype=object)),do_compression=True,long_field_names=True)
with (OUT/'raw_source_audit.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
