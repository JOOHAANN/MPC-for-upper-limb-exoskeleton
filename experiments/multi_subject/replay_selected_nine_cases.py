"""Re-run each subject's saved NMPC tuning without retraining BSE or MAP."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

from audit_subjects import ROOT

BASE=ROOT/'results/nine_case_75pct'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--subjects',nargs='+')
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    selected=json.loads((BASE/'selected_nmpc.json').read_text())
    mapping={
        '--guard-scale':'predictionErrorGuardScale',
        '--guard-floor':'predictionErrorGuardFloorNm',
        '--guard-cap':'predictionErrorGuardCapNm',
        '--under-weight':'backTorqueUnderWeight',
        '--effort-weight':'assistanceEffortWeight',
        '--headroom':'trackingHeadroomNm',
        '--safety-margin':'backTorqueMarginNm',
        '--horizon':'torqueConstraintHorizonSeconds',
        '--angle-weight':'angleWeight',
        '--velocity-weight':'velocityWeight',
        '--forecast':'inputForecast',
        '--forecast-time':'inputForecastTime',
        '--innovation-decay':'mapInnovationDecaySeconds',
    }
    for subject,candidate in selected.items():
        if args.subjects and subject not in args.subjects:
            continue
        report=json.loads((BASE/subject/'nmpc'/candidate/'evaluation.json').read_text())
        options=report['options']
        command=[sys.executable,str(Path(__file__).with_name('control_nine_cases.py')),
                 subject,'--candidate',candidate+'_replay','--workers',str(args.workers)]
        for flag,key in mapping.items():command.extend([flag,str(options[key])])
        if options.get('predictionErrorGuardRelative'):command.append('--relative-guard')
        if options.get('mapPhaseObserver'):command.append('--phase-observer')
        subprocess.run(command,cwd=ROOT,check=True)
        subprocess.run([sys.executable,str(Path(__file__).with_name('validate_nine_cases.py')),
                        '--subject',subject,'--candidate',candidate+'_replay'],cwd=ROOT,check=True)


if __name__=='__main__':main()
