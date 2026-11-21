"""Verify unchanged models and publish the single active comparison protocol."""
from pathlib import Path
import csv, hashlib, json
from datetime import datetime, timezone
import numpy as np
from scipy.io import loadmat

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
DEMO=ROOT/'Demo/NMPC'
OUT=DEMO/'output/pid_nmpc_528'
ACTIVE=OUT/'threshold_tracking'

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read_csv(path):
    with path.open(newline='',encoding='utf-8-sig') as f: return list(csv.DictReader(f))

def main():
    protected=json.loads((ROOT/'experiments/controller_comparison_528/original_file_hashes.json').read_text())
    frozen=json.loads((ROOT/'results/nine_case_75pct/release_manifest.json').read_text())
    files={**frozen['artifactSHA256'],**frozen['sourceCodeSHA256']}
    assert all(digest(ROOT/p)==v for p,v in protected.items())
    assert all(digest(ROOT/p)==v for p,v in files.items())
    log=(ROOT/'experiments/pid_damping_528/matlab_verification.log').read_text(errors='replace')
    assert 'PASS: one PID across 18 cases' in log
    assert 'PASS: clipped references, fixed tunings, shared 80-N-m caps, 200-N-m/s increase, instantaneous withdrawal' in log
    rows=read_csv(ACTIVE/'comparison_metrics_528.csv')
    summary=read_csv(ACTIVE/'comparison_summary_528.csv')
    assert len(rows)==54 and len(summary)==2
    result=loadmat(ACTIVE/'pid_nmpc_comparison_528.mat',simplify_cells=True)['results']
    assert result['mode']=='threshold_tracking'
    assert result['protocol']['torqueTarget']=='min(unassisted BSE torque, displayed muscle threshold)'
    assert result['algorithmSHA256']==digest(DEMO/'simulate_threshold_nmpc.m')
    audit=[]; pose=[]
    for group,g in enumerate(result['groups'],1):
        for source,pair in enumerate(g['runs'],1):
            original=loadmat(OUT/f'nmpc_group{group}_source{source}.mat',simplify_cells=True)
            assert np.max(np.abs(pair['baseline']-original['baseline']))<1e-7
            assert np.max(np.abs(pair['threshold']-original['threshold']))<1e-8
            target=np.minimum(pair['baseline'][:,3],pair['threshold'])
            assert np.max(np.abs(pair['nmpc']['target']-target))<1e-8
            angles=pair['nmpc']['actual'][:,0]
            pose.append(dict(Source=('Experimental activation','MAP activation')[source-1],Group=group,
                MinimumThetaRad=float(np.min(angles)),MaximumThetaRad=float(np.max(angles)),
                MaximumWholeSequenceTorqueExcessNm=float(np.max(pair['nmpc']['actual'][:,3]-pair['threshold'])),
                RecoveryAngleBelow1p8Rad=bool(np.max(angles)<1.8)))
            for method in ('pid','nmpc'):
                u=pair[method]['actual'][:,24]
                expected=np.minimum(pair[method]['command'][1:],u[:-1]+200*g['sequence']['dt'])
                assert np.max(np.abs(u[1:]-expected))<1e-8
                assert np.max(np.diff(u))/g['sequence']['dt']<=200+1e-6
            o=pair['options']
            assert o['instantRelease'] and not o['instantActuator'] and o['rateNmPerSec']==200
            assert o['backTorqueMarginNm']==.5 and o['predictionErrorGuardFloorNm']==0
            assert o['predictionErrorGuardScale']==0
    for source in ('Experimental activation','MAP activation'):
        for method in ('PID','NMPC'):
            selected=[r for r in rows if r['Source']==source and r['Controller']==method]
            assert len(selected)==9
            peak=max(float(r['AssistPeakNm']) for r in selected)
            rate=max(float(r['MaximumAssistRateNmPerSec']) for r in selected)
            assert peak<=80+1e-7

            sourceIndex=('Experimental activation','MAP activation').index(source)
            key='pid' if method=='PID' else 'nmpc'
            increase=max(float(np.max(np.diff(g['runs'][sourceIndex][key]['actual'][:,24]))/g['sequence']['dt']) for g in result['groups'])
            withdrawal=max(float(np.max(-np.diff(g['runs'][sourceIndex][key]['actual'][:,24]))) for g in result['groups'])
            assert increase<=200+1e-6
            audit.append(dict(Protocol='threshold_tracking',Source=source,Controller=method,
                SharedActuatorLimitNm=80,DeliveredPeakNm=peak,
                SharedTorqueIncreaseLimitNmPerSec=200,DeliveredMaximumIncreaseNmPerSec=increase,
                MaximumSingleStepWithdrawalNm=withdrawal,
                ActuatorResponse='increase limited to 200 N m/s; instantaneous withdrawal',DeliveredMaximumRateNmPerSec=rate))
    with (OUT/'assistance_limit_audit.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0])); w.writeheader(); w.writerows(audit)
    with (OUT/'pose_diagnostics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(pose[0])); w.writeheader(); w.writerows(pose)
    posePassed=all(r['RecoveryAngleBelow1p8Rad'] for r in pose)
    nmpcRows=[r for r in rows if r['Controller']=='NMPC']
    passed=sum(float(r['ExcessPeakNm'])<=1e-8 for r in nmpcRows)
    torquePassed=passed==18 and all(r['MaximumWholeSequenceTorqueExcessNm']<=1e-8 for r in pose)
    report=['# Subject-528 threshold-tracking PID/NMPC results','',
        'The active comparison uses one fixed NMPC tuning per activation source across all nine conditions. The PID has been retuned for damping with rise-rate-aware anti-windup; BSE/MAP models and NMPC settings are unchanged. RQ1 and RQ2 use the same closed-loop trajectories. There is no matched-effort constraint.','',
        'The NMPC torque target is exactly min(unassisted BSE torque, displayed threshold). The signed tracking error is penalized throughout actions and standing waits, including below-threshold valleys. Safety reserves apply to feasibility constraints only; they do not shift the reference. Figures include this upper-clipped target explicitly.','',
        'Both controllers share a 0-80-N-m actuator with torque increase limited to 200 N m/s and instantaneous withdrawal. The rising ramp is preserved inside integration substeps; the falling ramp is removed. Controller updates remain 0.02 s. Muscle activation and joint dynamics remain unchanged. The total prediction reserve is exactly 0.5 N m; additional prediction-error guards are disabled.','',
        'The actuator/reference contract is checked separately from torque-threshold and recovery-angle performance. Reduced reserves do not guarantee zero actual exceedance. See the measured counts and pose diagnostics below; no curves are clipped.','',
        '| Input | NMPC cases within threshold | Peak headroom range (N m) | PID/NMPC total excess area (N m s) | PID/NMPC mean angle deviation (deg) |',
        '| --- | ---: | ---: | ---: | ---: |']
    for r in summary:
        report.append(f"| {r['Source']} | {sum(float(x['ExcessPeakNm'])<=1e-8 for x in nmpcRows if x['Source']==r['Source'])}/9 | {float(r['NMPCMinimumPeakHeadroomNm']):.3f} to {float(r['NMPCMaximumPeakHeadroomNm']):.3f} | {float(r['PIDTotalExcessAreaNms']):.3f} / {float(r['NMPCTotalExcessAreaNms']):.3f} | {float(r['PIDMeanThetaRMSEDeg']):.3f} / {float(r['NMPCMeanThetaRMSEDeg']):.3f} |")
    report+=['',
        'Peak headroom is displayed threshold minus the maximum assisted back torque during each action. Exceedance counts are computed from the new asymmetric-actuator traces, not inherited from earlier runs. This does not imply a perfectly flat torque plateau: mean demand headroom is reported separately in the CSV.','',
        'No extra slow withdrawal is used for either activation source. Both controllers limit torque increase identically and immediately apply requested torque decreases. Standing waits retain the same NMPC clipped-reference objective. There is no separate fixed recovery-torque target or future measured-reference input.','',
        'The total internal feasibility reserve is 0.5 N m, without any additive uncertainty reserve. Both sources optimize 8 control blocks over 0.4 s. Prediction and actual output share a 0.01-s grid and 0.001-s inner integration. Experimental prediction covers 40 points (0.4 s). MAP adds a 0.8-s terminal zero-command recovery preview, for 120 points and 1.2 s total. Terminal zero commands use instantaneous withdrawal. Settings are fixed across all nine conditions within each source.','',
        'The active tunings use angle/velocity penalties inside the NMPC objective, without post-QP angle filtering. The optimized control sequence is retained and shifted between updates. No joint states are clipped and no extra restoring torques are added.','',
        'This run retunes only the reactive PID and adds rise-rate-aware conditional integration. BSE parameters, input data, NMPC settings, clipped-reference definition and actuator assumptions remain fixed. See pid_damping_verification.json for the before/after ripple diagnostics. Instantaneous withdrawal is an idealized assumption, not a claim of measured hardware response.','',
        '| Input | Load (kg) | Speed code | Threshold (N m) | NMPC back peak (N m) | Headroom (N m) | NMPC assistance peak (N m) |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    selected=sorted([r for r in rows if r['Controller']=='NMPC'],key=lambda r:(int(float(r['SourceIndex'])),float(r['LoadKg']),float(r['SpeedCode'])))
    for r in selected:
        report.append(f"| {r['Source']} | {float(r['LoadKg']):g} | {float(r['SpeedCode']):g} | {float(r['ThresholdNm']):.3f} | {float(r['BackPeakNm']):.3f} | {float(r['PeakHeadroomNm']):.3f} | {float(r['AssistPeakNm']):.3f} |")
    report+=['',
        'These are retrospective checks on the nine calibrated recorded scenarios, not independent validation or a hardware safety guarantee. Smaller angle deviation is not, by itself, a measure of synchronization, hardware latency or comfort. No claim of improved timing is made unless the separate timing results support it.','',
        'Active figures and raw data: `threshold_tracking/`. Earlier `actuator_matched/` and `effort_matched/` folders are historical and are not loaded by the current entry script. The equal-effort protocol is retired.','',
        f'Preservation checks passed: {len(protected)} original Demo source/data files and {len(files)} frozen-release files are unchanged. MATLAB verified shared caps, rate-limited increase and immediate withdrawal, the 0.5-N-m reserve, fixed source-specific controllers, torque accounting and causal forecasts. Threshold and pose performance are reported separately.']
    (OUT/'comparison_report.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    current=[DEMO/name for name in ['offline_PID_NMPC_comparison_528.m','simulate_threshold_nmpc.m',
        'simulate_threshold_pid.m','pid_nmpc_metrics.m','plot_pid_nmpc_comparison.m','README_PID_NMPC.md',
        'data/threshold_tracking_settings_528.mat','data/threshold_tracking_controller_settings.csv',
        'instant_bse_actuator.m','private/native/bse_instant_preview.c']]
    current.append(ROOT/'tests/test_pid_nmpc_comparison.m')
    for p in current:
        if p.suffix=='.m': assert all(ord(c)<128 for c in p.read_text(encoding='utf-8-sig'))
    verification=dict(createdAtUTC=datetime.now(timezone.utc).isoformat(),subject='240528Pilot',
        activeProtocol='threshold_tracking',effortMatched=False,conditionsPerSource=9,
        torqueTarget='min(unassisted BSE torque, displayed muscle threshold)',
        fullCycleSignedTorqueTracking=True,
        activationSources=2,passedNMPCActionChecks=passed,sameTracesForRQ1AndRQ2=True,
        samePIDAllCases=True,oneNMPCPerActivationSource=True,actuatorCapNm=80,
        torqueIncreaseRateNmPerSec=200,torqueDecreaseRateNmPerSec=None,instantaneousWithdrawal=True,totalTorqueReserveNm=.5,controlPeriodSeconds=.02,thresholdFraction=.75,
        originalDemoFilesUnchanged=len(protected),frozenReleaseFilesUnchanged=len(files),
        actuatorReferenceContractPassed=True,torqueThresholdChecksPassed=torquePassed,poseSanityPassed=posePassed,
        matlabVerificationPassed=(posePassed and torquePassed),overallStatus=('passed' if posePassed and torquePassed else 'contract_passed_performance_unresolved'),
        poseDiagnostics=pose,baselineInputsAndReferencesUnchanged=True,PIDGainsRetunedForDamping=True,NMPCTracesUnchanged=True,
        newFileSHA256={str(p.relative_to(ROOT)):digest(p) for p in current},summary=summary)
    (OUT/'comparison_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps({k:verification[k] for k in ['activeProtocol','effortMatched','passedNMPCActionChecks',
        'originalDemoFilesUnchanged','frozenReleaseFilesUnchanged','actuatorCapNm']},indent=2))

if __name__=='__main__': main()
