"""Search body geometry and neural policy jointly outside a local fit basin.

All scores use the same nine immutable conditions, references and MAP profiles.
This writes a candidate seed only, never a released or accepted model.
"""

from concurrent.futures import ThreadPoolExecutor
import json
import time

import numpy as np
from scipy.io import loadmat,savemat
from scipy.optimize import differential_evolution
from scipy.stats import qmc

import recalibrate_bse_pearson as calibrate


def main():
    subject='240730Exp_A01'; source=0
    package,models,groups,_,_=calibrate.load_subject(subject)
    folder=calibrate.OUT/subject/'source1'
    current=np.asarray(loadmat(folder/'checkpoint.mat',simplify_cells=True)['vector'])
    current=np.r_[current,np.zeros(48-len(current))]
    low=np.r_[calibrate.fit.LOW,calibrate.fit.NLOW,[-12,-12,-4,-6,-8,-.75]]
    high=np.r_[calibrate.fit.HIGH,calibrate.fit.NHIGH,[12,12,4,6,8,1.25]]
    low[2]=.65;low[38]=0;high[40]=1.5;high[37]=80
    span=high-low
    selected=np.r_[[0,5,8,12,13,18,19],np.arange(32,48)]
    center=np.clip((current-low)/span,1e-5,1-1e-5)
    records=[calibrate.fit.record(g['sequence'],models[source]['generator'],3) for g in groups]
    calibrate.extra_library(16)
    best=[np.inf,current.copy(),None]
    started=time.monotonic(); generation=[0]

    def objective(z):
        point=center.copy();point[selected]=z
        vector=low+span*point
        rows,_=calibrate.evaluate(vector[:32],vector[32:],package['controller'],records,source,6)
        stats=calibrate.summary(rows)
        value=(30*max(0,.82-stats['meanTorqueR']))**2+(40*max(0,.725-stats['meanThetaR']))**2
        value+=(12*max(0,stats['maximumAbdomenBackRatio']-.35))**2
        value+=(10*max(0,stats['maximumTheta']-1.95)/.2)**2
        value+=(10*max(0,-.10-stats['minimumTheta'])/.2)**2
        value+=(10*max(0,stats['boundaryPeakNm']-25)/25)**2
        value+=1000*stats['failed']+2*(1-stats['meanTorqueR'])+.6*(1-stats['meanThetaR'])
        value+=.2*(stats['meanTorqueRMSE']/80)**2
        if value<best[0]:best[:]=[float(value),vector.copy(),stats]
        return value

    def callback(x,convergence):
        generation[0]+=1
        savemat(folder/'global_seed.mat',dict(vector=best[1],loss=best[0]))
        print('GLOBAL',generation[0],'best',round(best[0],4),
            'rT/rQ',round(best[2]['meanTorqueR'],4),round(best[2]['meanThetaR'],4),
            'seconds',round(time.monotonic()-started),flush=True)
        return bool(best[2]['meanTorqueR']>.825 and best[2]['meanThetaR']>.73 and
            calibrate.quality_passed(best[2]))

    count=6*len(selected)
    rng=np.random.default_rng(73021)
    population=np.clip(center[selected]+rng.normal(0,.15,(count,len(selected))),1e-6,1-1e-6)
    population[:count//4]=qmc.LatinHypercube(len(selected),seed=730).random(count//4)
    population[0]=center[selected]
    # Explicitly include both gravity-dominated and load-dominated body seeds.
    for k,(body,lever,load) in enumerate([(15,.20,.85),(25,.25,.85),(30,.20,.70),(45,.35,.35)],1):
        vector=current.copy();vector[18]=body;vector[19]=lever;vector[12]=load
        population[k]=np.clip(((vector-low)/span)[selected],1e-6,1-1e-6)
    with ThreadPoolExecutor(max_workers=8) as pool:
        differential_evolution(objective,[(0.,1.)]*len(selected),init=population,
            workers=pool.map,updating='deferred',callback=callback,polish=False,
            maxiter=45,tol=1e-5,mutation=(.3,.9),recombination=.85,seed=73021)
    savemat(folder/'global_seed.mat',dict(vector=best[1],loss=best[0]))
    (folder/'global_search.json').write_text(json.dumps(dict(loss=best[0],metrics=best[2],
        vector=best[1].tolist(),selectedIndices=(selected+1).tolist(),referenceChanged=False),indent=2)+'\n')


if __name__=='__main__':
    main()
