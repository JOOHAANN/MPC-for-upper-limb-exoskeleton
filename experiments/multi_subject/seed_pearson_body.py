"""Screen distinct body/load-dominated initial conditions on all nine cases."""
from concurrent.futures import ThreadPoolExecutor
from itertools import product
import json
import numpy as np
from scipy.io import savemat
import recalibrate_bse_pearson as c


def main():
    subject='240730Exp_A01';source=0
    package,models,groups,_,_=c.load_subject(subject)
    records=[c.fit.record(g['sequence'],models[0]['generator'],3) for g in groups]
    c.extra_library(16)
    original=models[0]['parameters']
    choices=list(product([15,25,40],[.20,.30],[.60,.85],[.75,1.],[.3,.8,1.5],[.9,1.4],[0,.5]))
    def run(choice):
        mass,lever,load,rho,lead,amplitude,tonic=choice
        p=original.copy();p[18]=mass;p[19]=lever;p[12]=load;p[13]=1.1
        p[7]=.02;p[8]=200;p[9]=1.;p[14:18]=[.4,.02,.05,.05]
        r=np.r_[[amplitude,.7,.05,160,35,60,rho,0,lead,tonic],np.zeros(5),.15]
        rows,_=c.evaluate(p,r,package['controller'],records,source,6)
        s=c.summary(rows)
        cost=(30*max(0,.82-s['meanTorqueR']))**2+(40*max(0,.725-s['meanThetaR']))**2
        cost+=(12*max(0,s['maximumAbdomenBackRatio']-.35))**2
        cost+=(10*max(0,s['maximumTheta']-1.95)/.2)**2
        cost+=(10*max(0,-.1-s['minimumTheta'])/.2)**2
        cost+=(10*max(0,s['boundaryPeakNm']-25)/25)**2+1000*s['failed']
        return cost,np.r_[p,r],s,choice
    with ThreadPoolExecutor(max_workers=6) as pool:
        results=list(pool.map(run,choices))
    results.sort(key=lambda item:item[0])
    folder=c.OUT/subject/'source1'
    savemat(folder/'body_seed.mat',dict(vector=results[0][1]))
    (folder/'body_seed_audit.json').write_text(json.dumps([
        dict(cost=cost,metrics=metrics,choice=choice,vector=vector.tolist())
        for cost,vector,metrics,choice in results[:20]],indent=2)+'\n')
    print(json.dumps(dict(candidates=len(choices),cost=results[0][0],metrics=results[0][2],choice=results[0][3]),indent=2))


if __name__=='__main__':
    main()
