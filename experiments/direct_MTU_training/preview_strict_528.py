"""Inspect raw strict-input candidates before publishing the MATLAB replay."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import strict_fit_528 as f

subject=f.loadmat(f.HERE/'strict_input_528.mat',simplify_cells=True)['subject']
paths=['strict_global_528_1_seed1.json','strict_528_2.json']
reports=[json.loads((f.HERE/p).read_text()) for p in paths]
rr=[f.records(subject,s,.01) for s in range(2)]
outs=[[f.rollout(np.ascontiguousarray(r['parameters']),r['abdomenScale'],q) for q in rr[s]] for s,r in enumerate(reports)]
fig,axes=plt.subplots(4,3,figsize=(15,13));colors=['#0069b4','#8c3896']
q=f.cells(subject['sequences'])[0]
for k,seg in enumerate(f.cells(q['segments'])):
    a=int(seg['startIndex'])-1;b=int(seg['endIndex']);ix=np.arange(max(0,a-50),min(len(q['time']),b+100));t=q['time'][ix]-q['time'][a]
    axes[0,k].plot(t,q['referenceTorque'][ix],'r--',lw=1,label='CEINMS reference')
    axes[1,k].plot(t,q['referenceTheta'][ix],'r--',lw=1,label='Original IK reference')
    for s,out in enumerate(outs):
        x=out[0][ix]
        axes[0,k].plot(t,x[:,2],color=colors[s],label=['Experimental','MAP'][s])
        axes[1,k].plot(t,x[:,0],color=colors[s])
        axes[2,k].plot(t,x[:,3],color=colors[s],label=['Exp. back','MAP back'][s])
        axes[2,k].plot(t,x[:,4],'--',color=colors[s],label=['Exp. abdomen','MAP abdomen'][s])
        axes[3,k].plot(t,x[:,21],color=colors[s],label=['Exp. standing PD','MAP standing PD'][s])
    axes[3,k].plot(t,outs[0][0][ix,25],'k--',label='Boundary torque: zero')
    axes[0,k].set_title(f"{seg['loadKg']:g} kg | code {seg['speedCode']:g}")
    for row in range(4):axes[row,k].grid(alpha=.25)
    axes[3,k].set_xlabel('Time from recorded action onset (s)')
for row,ylabel in enumerate(['Net torque (N m)','Trunk angle (rad)','Muscle torque (N m)','Standing PD request (N m)']):
    axes[row,0].set_ylabel(ylabel);axes[row,0].legend(fontsize=8)
fig.suptitle('Subject 528 | Strict experimental/independently trained MAP abdomen | Candidate preview\nNo movement PD, assistance or angle boundaries')
fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(f.HERE/'strict_preview.png',dpi=130)
print('Saved strict_preview.png')
