"""Rebuild selected references and activations from unchanged source files."""
from __future__ import annotations

import json
import numpy as np
from scipy.io import loadmat

from audit_subjects import ROOT, DATA, read_sto
from control_nine_cases import BASE, digest, many
from prepare_exp_three_loads import angle_reference
from recompute_exp_activation import selected_signal, MODEL_CHANNELS
from prepare_subjects import map_predict
from quality_gate import gaussian_matlab


def main():
    rows=[]
    for subject in ('240725Exp_A01','240730Exp_A01'):
        folder=BASE/subject
        data=loadmat(folder/'prepared.mat',simplify_cells=True)
        gate=json.loads((folder/'quality_gate.json').read_text())
        for cycle,audit in zip(many(data['cycles']),gate['trials']):
            trial=DATA/subject/cycle['trialPath']
            for name,expected in audit['sourceFiles'].items():
                assert digest(trial/name)==expected, 'Source changed: '+str(trial/name)
            origin=float(read_sto(trial/'emgFilt.sto')[1][0,0])
            epoch=origin+cycle['sourceStartSeconds']+.01*np.arange(len(cycle['q']))
            angle,_=angle_reference(trial,epoch)
            error=float(np.max(np.abs(angle-cycle['q'])))
            assert error<2e-5, (subject,cycle['id'],error)
            selected=dict(startEpochSeconds=float(epoch[0]-1e-6),
                          endEpochSeconds=float(epoch[-1]+1e-6),
                          source='CEINMS_activation' if audit['source']=='model-input EMG' else 'filtered')
            row=dict(subject=subject.removesuffix('_A01'),trial=cycle['trialPath'],
                     selected=selected,originEpochSeconds=origin)
            signal,_=selected_signal(row,{name:.1 for name in MODEL_CHANNELS},
                    additional_smoothing_ms=audit['additionalSubjectWideGaussianWindowMs'])
            back=signal['activation'][:,:6]@data['mapModel']['lumpWeights']
            activation_error=float(np.max(np.abs(back-cycle['back'])))
            assert activation_error<1e-10
            generated=map_predict(data['mapModel'],cycle)
            if audit['additionalSubjectWideGaussianWindowMs']:
                generated=gaussian_matlab(generated,31)
            map_error=float(np.max(np.abs(generated-cycle['mappedBack'])))
            assert map_error<1e-10
            valid=np.asarray(cycle['torqueReferenceValid'],dtype=bool)
            assert valid.any() and np.isnan(cycle['torque'][~valid]).all()
            rows.append(dict(subject=subject,recordId=cycle['id'],unchangedSourceFiles=len(audit['sourceFiles']),
                originalScaleReferenceMaximumDifferenceRad=error,
                experimentalActivationMaximumDifference=activation_error,
                trainedMAPActivationMaximumDifference=map_error,
                availableTorqueFraction=float(valid.mean())))
    path=BASE/'input_audit.json'
    path.write_text(json.dumps(dict(passed=True,rows=rows,noSelectedBoutRescaling=True,
                        allConditionsHaveObservedTorqueSupport=True),indent=2)+'\n')
    print(json.dumps(dict(passed=True,conditions=len(rows),
                    maximumReferenceDifferenceRad=max(r['originalScaleReferenceMaximumDifferenceRad'] for r in rows)),indent=2))


if __name__=='__main__':main()
