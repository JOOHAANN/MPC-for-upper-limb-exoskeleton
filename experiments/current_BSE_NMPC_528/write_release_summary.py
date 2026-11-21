"""Summarize the saved MATLAB replay without modifying controller traces."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FOLDER = ROOT / 'Demo/NMPC'
OUT = FOLDER / 'output/current_BSE_528'
with (OUT / 'comparison_metrics_528.csv').open(newline='') as file:
    rows = list(csv.DictReader(file))
settings = loadmat(FOLDER / 'data/current_controller_528.mat', simplify_cells=True)['settings']
report = {'subject': '240528Pilot', 'conditionsPerSource': 9,
          'scope': 'Offline calibration-set replay; not a hardware or unseen-condition validation.',
          'modelSHA256': hashlib.sha256((ROOT / 'Demo/BSE/Direct_Activation_No_Limits/data/subject_528.mat').read_bytes()).hexdigest(),
          'sources': {}}
for source, name in enumerate(['Experimental activation', 'MAP activation'], 1):
    group_rows = [r for r in rows if int(r['SourceIndex']) == source]
    baseline = [r for r in group_rows if r['Controller'] == 'Unassisted']
    nmpc = [r for r in group_rows if r['Controller'] == 'NMPC']
    pid = [r for r in group_rows if r['Controller'] == 'PID']
    total = sum(float(r['ExcessAreaNms']) for r in baseline)
    flags = []
    corrections = []
    for group in range(1, 4):
        d = loadmat(OUT / f'group{group}_source{source}.mat', simplify_cells=True)
        n = d['pair']['nmpc']
        flags.extend(n['solverExitFlag'][np.isfinite(n['solverExitFlag'])].tolist())
        # Columns 20-22 are zero only during movement, as verified by MATLAB.
        moving = n['actual'][:, 19] == 0
        corrections.extend((n['actual'][moving, 17] - d['inputs'][moving, 0]).tolist())
    worst = max(nmpc, key=lambda r: float(r['ExcessPeakNm']))
    report['sources'][name] = dict(
        maximumNMPCExcessNm=float(worst['ExcessPeakNm']),
        worstCase=dict(loadKg=float(worst['LoadKg']), speedCode=int(worst['SpeedCode'])),
        conditionsAboveThreshold=sum(float(r['ExcessPeakNm']) > 1e-6 for r in nmpc),
        conditionsAboveThresholdByMoreThanHalfNm=sum(float(r['ExcessPeakNm']) > .5 for r in nmpc),
        nmpcExcessAreaReductionPercent=100*(1-sum(float(r['ExcessAreaNms']) for r in nmpc)/total),
        pidExcessAreaReductionPercent=100*(1-sum(float(r['ExcessAreaNms']) for r in pid)/total),
        meanNMPCThetaDeviationDeg=float(np.mean([float(r['ThetaDeviationRMSEDeg']) for r in nmpc])),
        maximumNMPCActuatorTorqueNm=max(float(r['AssistPeakNm']) for r in nmpc),
        maximumNMPCPositiveRiseRateNmPerSec=max(float(r['MaximumAssistRiseRateNmPerSec']) for r in nmpc),
        maximumBoundaryTorqueNm=max(float(r['BoundaryPeakNm']) for r in group_rows),
        unsuccessfulUpdateCount=sum(f <= 0 for f in flags),
        updatesAbove20ms=sum(int(r['UpdatesAbove20ms']) for r in nmpc),
        totalUpdates=sum(int(r['TotalUpdates']) for r in nmpc),
        maximumRecordedSolveMs=max(float(r['SolveMaxMs']) for r in nmpc),
        backCommandCorrectionRange=[min(corrections), max(corrections)])
fields = ['Source', 'OuterSampleIntervalSeconds', 'InternalStepSeconds',
          'ControllerUpdateSeconds', 'PredictionPoints', 'PredictionStepSeconds',
          'PredictionHorizonSeconds', 'ControlledPoints', 'ControlBlocks',
          'ControlWindowSeconds', 'ZeroCommandRecoveryTailSeconds',
          'AssistanceMinimumNm', 'AssistanceMaximumNm', 'PositiveRiseLimitNmPerSec',
          'WithdrawalRateLimitNmPerSec', 'ThresholdFraction', 'FeasibilityReserveNm',
          'MovementPD_P', 'MovementPD_D', 'PDOnlyWithAppliedAssistance']
with (FOLDER / 'data/current_controller_settings_528.csv').open('w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(fields)
    for source, name in enumerate(['Experimental activation', 'MAP activation']):
        o = settings['controllers'][source]
        writer.writerow([name, .01, .001, o['controlPeriod'], o['points'], o['predictionDt'],
                         o['points']*o['predictionDt'], o['controlPoints'], o['blocks'],
                         o['controlPoints']*o['predictionDt'], (o['points']-o['controlPoints'])*o['predictionDt'],
                         0, o['capNm'], o['rateNmPerSec'], 'Inf', .75, .5, 5, 3, True])
(OUT / 'release_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2))
