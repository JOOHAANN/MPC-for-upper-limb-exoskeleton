"""Audit and plot the delivered MATLAB PID traces; never smooth output curves."""
from pathlib import Path
import csv, json
import numpy as np
from scipy.io import loadmat
from scipy.signal import butter, sosfiltfilt, find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
OUT=ROOT/'Demo/NMPC/output/pid_nmpc_528'
ACTIVE=OUT/'threshold_tracking'
FILTER=butter(3,3,fs=100,btype='highpass',output='sos')

def diagnostics(out,idx):
    u=out[idx,24]; back=out[idx,3]
    return dict(AssistanceRippleRMSNm=float(np.sqrt(np.mean(sosfiltfilt(FILTER,u)**2))),
        BackRippleRMSNm=float(np.sqrt(np.mean(sosfiltfilt(FILTER,back)**2))),
        AssistanceTotalVariationNm=float(np.sum(np.abs(np.diff(u)))),
        AssistanceReversals=int(len(find_peaks(u,prominence=2)[0])+len(find_peaks(-u,prominence=2)[0])))

def main():
    old=loadmat(HERE/'before_comparison_528.mat',simplify_cells=True)['results']
    new=loadmat(ACTIVE/'pid_nmpc_comparison_528.mat',simplify_cells=True)['results']
    selection=json.loads((HERE/'selection.json').read_text())
    rows=[]; cases=[]; maximum_parity=0.
    # Check the native development rollout against the final MATLAB run.
    import tune_damping as tuning
    prepared=tuning.previous.prepare()
    for group,g in enumerate(new['groups']):
        q=g['sequence']
        for source,pair in enumerate(g['runs']):
            previous=old['groups'][group]['runs'][source]
            assert np.array_equal(pair['nmpc']['actual'],previous['nmpc']['actual']), 'NMPC changed'
            assert np.array_equal(pair['baseline'],previous['baseline']), 'Baseline changed'
            assert np.array_equal(pair['threshold'],previous['threshold']), 'Threshold changed'
            assert pair['pid']['pid']==new['pid']
            c=prepared[group*2+source]
            c['a']=np.ascontiguousarray(pair['nmpc']['inputs'])
            c['threshold']=np.ascontiguousarray(pair['threshold'])
            native=tuning.previous.roll(c,np.array(selection['gains']))
            parity=float(np.max(np.abs(native-pair['pid']['actual'])))
            maximum_parity=max(maximum_parity,parity)
            assert parity<.001, f'Development/MATLAB discrepancy {parity}'
            for stage,seg in enumerate(q['segments']):
                idx=slice(int(seg['startIndex'])-1,int(seg['endIndex']))
                for version,p in [('Previous',previous),('Damped',pair)]:
                    row=dict(Source=('Experimental activation','MAP activation')[source],Group=group+1,
                        LoadKg=seg['loadKg'],SpeedCode=seg['speedCode'],PID=version)
                    row.update(diagnostics(p['pid']['actual'],idx))
                    row['ExcessPeakNm']=float(max(0,np.max(p['pid']['actual'][idx,3]-pair['threshold'][idx])))
                    row['ExcessAreaNms']=float(np.sum(np.maximum(p['pid']['actual'][idx,3]-pair['threshold'][idx],0))*q['dt'])
                    rows.append(row)
                cases.append((source,q,pair,previous,seg,idx))
    with (OUT/'pid_damping_metrics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary={}
    for source in ('Experimental activation','MAP activation'):
        summary[source]={}
        for version in ('Previous','Damped'):
            selected=[r for r in rows if r['Source']==source and r['PID']==version]
            summary[source][version]={key:float(np.mean([r[key] for r in selected])) for key in selected[0]
                if key not in ('Source','Group','LoadKg','SpeedCode','PID')}
    before=[r for r in rows if r['PID']=='Previous']; after=[r for r in rows if r['PID']=='Damped']
    assert all(b['AssistanceRippleRMSNm']<a['AssistanceRippleRMSNm'] for a,b in zip(before,after))
    assert all(b['AssistanceTotalVariationNm']<a['AssistanceTotalVariationNm'] for a,b in zip(before,after))
    report=dict(gains=new['pid'],samePIDAll18Cases=True,NMPCTracesUnchanged=True,
        baselineAndThresholdUnchanged=True,maximumNativeMATLABDiscrepancy=maximum_parity,
        rippleReducedCases=18,totalVariationReducedCases=18,summary=summary,
        diagnostics='Third-order zero-phase high-pass at 3 Hz, applied only to diagnostic metrics. Displayed traces are unfiltered. Reversals use 2-N-m prominence. These are finite-trial ripple diagnostics, not a formal stability proof.')
    (OUT/'pid_damping_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    plt.rcParams.update({'font.family':'serif','font.size':9,'axes.grid':True,'grid.alpha':.2})
    for source in range(2):
        selected=[c for c in cases if c[0]==source]
        selected.sort(key=lambda c:(c[4]['loadKg'],c[4]['speedCode']))
        fig,axes=plt.subplots(3,3,figsize=(14,8),layout='constrained')
        for ax,(_,q,p,previous,seg,idx) in zip(axes.flat,selected):
            time=q['time'][idx]-q['time'][idx.start]
            ax.plot(time,previous['pid']['actual'][idx,24],color='.6',lw=1,label='Previous PID')
            ax.plot(time,p['pid']['actual'][idx,24],color='#df6e18',lw=1.3,label='Damped PID')
            ax.set_title(f"{seg['loadKg']:g} kg / code {seg['speedCode']:g}")
            ax.set(xlabel='Time within action (s)',ylabel='Assistance torque (N m)')
        axes.flat[0].legend()
        fig.suptitle('PID damping comparison | '+('Experimental activation','MAP activation')[source])
        stem=ACTIVE/('PID_damping_all_nine_'+('experimental','MAP')[source])
        fig.savefig(str(stem)+'.png',dpi=250); fig.savefig(str(stem)+'.pdf'); plt.close(fig)
    # Show the former worst-ripple action at full sample resolution.
    worst=max(before,key=lambda r:r['AssistanceRippleRMSNm'])
    source=('Experimental activation','MAP activation').index(worst['Source'])
    c=next(c for c in cases if c[0]==source and c[4]['loadKg']==worst['LoadKg'] and c[4]['speedCode']==worst['SpeedCode'])
    _,q,p,previous,seg,idx=c
    time=q['time'][idx]-q['time'][idx.start]
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True,layout='constrained')
    for ax,column in zip(axes,(3,24)):
        ax.plot(time,previous['pid']['actual'][idx,column],color='.6',lw=1,label='Previous PID')
        ax.plot(time,p['pid']['actual'][idx,column],color='#df6e18',lw=1.5,label='Damped PID')
        if column==3:
            ax.plot(time,p['baseline'][idx,column],'r--',lw=.8,label='Unassisted')
            ax.plot(time,p['threshold'][idx],'k--',lw=1,label='Muscle threshold')
        else: ax.axhline(80,color='.3',ls='--',lw=1,label='Shared actuator limit')
        ax.set_ylabel(('Back-muscle' if column==3 else 'Assistance')+' torque (N m)')
        ax.legend(ncol=4,fontsize=8)
    axes[-1].set_xlabel('Time within action (s)')
    fig.suptitle(f"PID damping | {worst['Source']} | {seg['loadKg']:g} kg / code {seg['speedCode']:g}")
    fig.savefig(ACTIVE/'PID_damping_example.png',dpi=250)
    fig.savefig(ACTIVE/'PID_damping_example.pdf')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
