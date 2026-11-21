"""Verify original-file preservation and publish an English result report."""
from pathlib import Path
import csv, hashlib, json
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
DEMO=ROOT/'Demo/NMPC'
OUT=DEMO/'output/pid_nmpc_528'

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    protected=json.loads((HERE/'original_file_hashes.json').read_text())
    assert all(digest(ROOT/name)==value for name,value in protected.items())
    frozen=json.loads((ROOT/'results/nine_case_75pct/release_manifest.json').read_text())
    frozen_files={**frozen['artifactSHA256'],**frozen['sourceCodeSHA256']}
    assert all(digest(ROOT/name)==value for name,value in frozen_files.items())
    log=(HERE/'verification_final.log').read_text(errors='replace')
    assert 'PASS: fixed PID' in log, 'Wait for the final MATLAB verification.'
    summary=list(csv.DictReader((OUT/'comparison_summary_528.csv').open()))
    assert len(summary)==4
    audit=[]
    for mode in ('actuator_matched','effort_matched'):
        rows=list(csv.DictReader((OUT/mode/'comparison_metrics_528.csv').open()))
        assert len(rows)==54
        for source in ('Experimental activation','MAP activation'):
            for method in ('PID','NMPC'):
                selected=[r for r in rows if r['Source']==source and r['Controller']==method]
                assert len(selected)==9
                cap=max(float(r['AssistPeakNm']) for r in selected)
                rate=max(float(r['MaximumAssistRateNmPerSec']) for r in selected)
                assert cap<=80+1e-7 and rate<=200+1e-5
                audit.append(dict(Protocol=mode,Source=source,Controller=method,
                    SharedActuatorLimitNm=80,DeliveredPeakNm=cap,
                    SharedRateLimitNmPerSec=200,DeliveredMaximumRateNmPerSec=rate))
    with (OUT/'assistance_limit_audit.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
    report=['# Subject-528 PID/NMPC comparison results','',
       'Two protocols are shown separately. PID gains, BSE parameters, inputs, thresholds and physical actuator constraints are fixed. All nine load/speed conditions and both activation sources are retained.','',
       '| Protocol | Input | PID effort | NMPC effort | Mismatch | PID/NMPC mean excess peak (N m) | PID/NMPC angle RMSE (deg) |',
       '| --- | --- | ---: | ---: | ---: | --- | --- |']
    for r in summary:
        report.append(f"| {r['Protocol']} | {r['Source']} | {float(r['PIDTotalEffortNm2s']):.2f} | {float(r['NMPCTotalEffortNm2s']):.2f} | {float(r['EffortMismatchPercent']):+.3f}% | {float(r['PIDMeanExcessPeakNm']):.3f} / {float(r['NMPCMeanExcessPeakNm']):.3f} | {float(r['PIDMeanThetaRMSEDeg']):.3f} / {float(r['NMPCMeanThetaRMSEDeg']):.3f} |")
    report += ['',
       'Effort is the integral of squared delivered assistance torque (N^2 m^2 s), not energy. Matching is across all nine conditions per source, not within every individual case. Angle RMSE is relative to the unassisted BSE trajectory.','',
       'The actuator-matched protocol retains the original NMPC options. It reduces threshold exceedance, but uses more effort and increases mean angle deviation. The effort-matched protocol changes the NMPC trade-off only; it retains the same PID and physical constraints. It reduces mean angle deviation but does not universally reduce torque exceedance.','',
       '| Effort-matched input | PID absolute demand lag (ms) | NMPC absolute demand lag (ms) | NMPC angle-RMSE improvements |',
       '| --- | ---: | ---: | ---: |']
    for r in summary:
        if r['Protocol']=='effort_matched':
            report.append(f"| {r['Source']} | {float(r['PIDMeanAbsoluteDemandLagMs']):.2f} | {float(r['NMPCMeanAbsoluteDemandLagMs']):.2f} | {r['NMPCThetaWins']}/9 |")
    report += ['',
       'The measured demand-lag proxy does not establish better temporal alignment. Smaller angle deviation alone must not be presented as reduced hardware latency or better comfort. This is a retrospective in-silico study using fixed calibrated models and protocol-known load/wait schedules.','',
       'Primary torque figures: `actuator_matched/RQ1_torque_all_nine_*.pdf`. Primary matched-effort motion figures: `effort_matched/RQ2_theta_all_nine_*.pdf` and `effort_matched/RQ2_theta_deviation_all_nine_*.pdf`. The exact original-controller angle comparison is also retained under `actuator_matched`.','',
       'Both protocols include continuous 5/10/15-kg examples, all raw traces, case-level metrics and actuator-limit audits. No traces are shifted, amplitude-scaled or state-clipped after simulation.','',
       f'Verification passed: {len(protected)} pre-existing Demo source/data files and {len(frozen_files)} original frozen-release files are unchanged. Fixed gains, actuator cap/rate, control period, thresholds, effort matching, accounting and PID causality passed the MATLAB test.']
    (OUT/'comparison_report.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    new_code=[p for p in DEMO.glob('*.m') if str(p.relative_to(ROOT)) not in protected]
    files=new_code+[DEMO/'data/comparison_settings_528.mat',DEMO/'README_PID_NMPC.md',
                    ROOT/'tests/test_pid_nmpc_comparison.m']
    for p in new_code:
        assert not any('\u4e00'<=c<='\u9fff' for c in p.read_text(encoding='utf-8-sig'))
    manifest=dict(createdAtUTC=datetime.now(timezone.utc).isoformat(),subject='240528Pilot',
        protocols=['actuator_matched','effort_matched'],conditionsPerSource=9,activationSources=2,
        originalDemoFilesUnchanged=len(protected),frozenReleaseFilesUnchanged=len(frozen_files),
        originalModelAndDataUnchanged=True,independentPhysiologicalValidation=False,
        samePIDAllCasesAndProtocols=True,actuatorCapNm=80,actuatorRateNmPerSec=200,
        controllerPeriodSeconds=.02,thresholdFraction=.75,
        effortMatchScope='total over all nine conditions per activation source',
        originalFileSHA256=protected,newFileSHA256={str(p.relative_to(ROOT)):digest(p) for p in files},
        summary=summary,matlabVerificationPassed=True)
    (OUT/'comparison_verification.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ('originalDemoFilesUnchanged','frozenReleaseFilesUnchanged',
         'samePIDAllCasesAndProtocols','actuatorCapNm','matlabVerificationPassed')},indent=2))

if __name__=='__main__':main()
