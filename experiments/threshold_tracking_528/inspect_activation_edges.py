"""Inspect causal forecast scenarios with the exact delivered Hill kernel."""
from pathlib import Path
import ctypes,sys,numpy as np
from scipy.io import loadmat

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent.parent
sys.path.insert(0,str(ROOT/'experiments/controller_comparison_528'))
from tune_pid import inputs

def main():
    s=loadmat(ROOT/'Demo/BSE/data/subject_528.mat',simplify_cells=True)['subject']
    q=s['sequences'][2]; m=s['models'][0]; a=inputs(q,m,0); planned=inputs(q,m,1)
    cache=loadmat(ROOT/'Demo/NMPC/output/pid_nmpc_528/nmpc_group3_source1.mat',simplify_cells=True)
    options=cache['options']; pred=options['activationResidualPredictor']; lags=pred['lags'].astype(int)
    lib=ctypes.CDLL(str(ROOT/'experiments/controller_comparison_528/pid_tuning_core.dll'))
    arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    lib.hill_mpc_init.argtypes=[arr]*4+[ctypes.c_double]*2+[arr]*3
    lib.hill_mpc_advance.argtypes=[arr]*5+[ctypes.c_double]*4+[ctypes.c_int,arr]
    lib.hill_mpc_predict_map.argtypes=[arr]*8+[ctypes.c_int]*2+[ctypes.c_double,ctypes.c_int,ctypes.c_double,arr]
    p=np.ascontiguousarray(m['parameters']); r=np.ascontiguousarray(m['neuralParameters']); cfg=np.ascontiguousarray(np.r_[q['standingAngle'],s['balanceController']]); state=np.zeros(34); out=np.zeros(27)
    lib.hill_mpc_init(p,r,cfg,a[0],float(q['mass'][0]),float(q['waiting'][0]),q['initialState'],state,out)
    unc=loadmat(HERE/'activation_uncertainty.mat',simplify_cells=True)['uncertainty']; time=np.arange(1,17)*.025
    for current in range(1,1051):
        lib.hill_mpc_advance(p,r,cfg,state,a[current],float(q['mass'][current]),float(q['waiting'][current]),float(cache['run']['actual'][current,24]),.01,10,out)
        if current not in (1019,1029,1039): continue
        query=current+time/.01; idx=np.arange(len(a))
        activation=np.column_stack([np.interp(query,idx,planned[:,j]) for j in (0,1)])
        x=(a[np.maximum(current-lags,0)]-planned[np.maximum(current-lags,0)]).reshape(-1)/pred['scale']
        activation=np.clip(activation+(x@pred['coefficients']).reshape(16,2),0,1)
        mass=np.ascontiguousarray(np.interp(query,idx,q['mass'])); wait=np.ascontiguousarray(np.interp(query,idx,q['waiting']))
        print('TIME',current*.01)
        for strength in [.5,1,1.5]:
            for edge in [.05,.1,.15,.2,.25,.3]:
                change=strength*unc*np.tanh((time[:,None]-edge)/.025)*[1,-1]
                act=np.ascontiguousarray(np.clip(activation+change,0,1))
                forecast=np.zeros((16,27)); controls=np.zeros(4)
                flag=lib.hill_mpc_predict_map(p,r,cfg,state,controls,act,mass,wait,4,16,.025,10,200.,forecast)
                assert not flag
                print(strength,edge,round(forecast[:,3].max(),2))

if __name__=='__main__': main()
