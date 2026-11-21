from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'multi_subject'))
from audit_subjects import extract,DATA
fig,axes=plt.subplots(2,1,figsize=(14,7))
for ax,speed in zip(axes,[25,35]):
    q=extract(DATA/'240528Pilot'/f'dynamics_stoop_5kg_{speed}')
    ax.plot(q['t'],q['back'],label='Experimental torque-weighted back activation')
    peaks,_=find_peaks(q['back'],prominence=.02,distance=25)
    for i in peaks:ax.text(q['t'][i],q['back'][i],f"{q['t'][i]:.2f}",fontsize=8)
    ax.set_title(f'5 kg / code {speed} | time from common recorded support onset')
    ax.set_xlabel('Time (s)');ax.grid(alpha=.2);ax.legend()
fig.tight_layout();fig.savefig(HERE/'window_audit.png',dpi=150)
