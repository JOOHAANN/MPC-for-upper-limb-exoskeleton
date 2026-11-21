function rows = pid_nmpc_metrics(sequence,pair,group,source,sourceName)
%PID_NMPC_METRICS Distinguish torque tracking, motion change and timing.
% Angle error is relative to the unassisted BSE rollout, not the IK reference.
% Demand is the positive part of unassisted back torque minus the threshold.
% Positive correlation lag means delivered assistance follows that demand.
% Effort includes each action's allocated standing/release interval. These
% intervals partition the complete sequence, so their efforts sum exactly.
segments = sequence.segments;
count = numel(sequence.time);
cuts = [1,round((segments{1}.endIndex+segments{2}.startIndex)/2)+1, ...
    round((segments{2}.endIndex+segments{3}.startIndex)/2)+1,count+1];
controllers = {'Unassisted','PID','NMPC'};
outputs = {pair.baseline,pair.pid.actual,pair.nmpc.actual};
rows = struct([]);
for stage = 1:3
    segment = segments{stage};
    action = segment.startIndex:segment.endIndex;
    context = max(1,action(1)-40):min(count,action(end)+80);
    allocated = cuts(stage):cuts(stage+1)-1;
    threshold = pair.threshold(action(1));
    demand = max(0,pair.baseline(context,4)-threshold);
    for method = 1:3
        actual = outputs{method};
        excess = max(0,actual(action,4)-threshold);
        target=min(pair.baseline(action,4),threshold);
        trackingError=actual(action,4)-target;
        demandIndex=action(pair.baseline(action,4)>threshold);
        headroom=threshold-actual(demandIndex,4);
        error = actual(action,1)-pair.baseline(action,1);
        u = actual(context,25);
        [lag,correlation] = demand_lag(demand,u,sequence.dt);
        demandOnset = sustained_onset(demand,max(2,.05*max(demand)));
        assistOnset = sustained_onset(u,max(2,.05*max(u)));
        onsetOffset = (assistOnset-demandOnset)*sequence.dt*1000;
        row = struct('Subject',"240528Pilot",'Source',string(sourceName),'SourceIndex',source, ...
            'Group',group,'Stage',stage,'LoadKg',double(segment.loadKg), ...
            'SpeedCode',double(segment.speedCode),'Controller',string(controllers{method}), ...
            'ThresholdNm',threshold,'ActuatorLimitNm',pair.options.capNm, ...
            'ActuatorRateLimitNmPerSec',pair.options.rateNmPerSec, ...
            'ControlPeriodSeconds',pair.options.controlPeriod, ...
            'BackPeakNm',max(actual(action,4)),'ExcessPeakNm',max(excess), ...
            'ClippedTargetRMSEnm',sqrt(mean(trackingError.^2)), ...
            'ClippedTargetMAEnm',mean(abs(trackingError)), ...
            'PeakHeadroomNm',threshold-max(actual(action,4)), ...
            'DemandMeanHeadroomNm',mean(max(headroom,0)), ...
            'DemandWithin5NmPercent',100*mean(headroom>=0 & headroom<=5), ...
            'ExcessAreaNms',sum(excess)*sequence.dt, ...
            'ExcessDurationSeconds',sum(excess>.1)*sequence.dt, ...
            'ThetaDeviationRMSEDeg',rad2deg(sqrt(mean(error.^2))), ...
            'ThetaDeviationMaximumDeg',rad2deg(max(abs(error))), ...
            'ThetaRecoveryRMSEDeg',rad2deg(sqrt(mean((actual(context,1)-pair.baseline(context,1)).^2))), ...
            'ControlEffortNm2s',sum(actual(allocated,25).^2)*sequence.dt, ...
            'AssistPeakNm',max(actual(allocated,25)), ...
            'MaximumAssistRateNmPerSec',max(abs(diff(actual(allocated,25))))/sequence.dt, ...
            'DemandLagMs',lag,'DemandCorrelation',correlation,'OnsetOffsetMs',onsetOffset, ...
            'AbdomenPeakNm',max(actual(action,5)),'BoundaryPeakNm',max(abs(actual(context,26))));
        rows = [rows;row]; %#ok<AGROW>
    end
end
end

function [lagMs,best] = demand_lag(demand,assistance,dt)
lagMs = NaN; best = NaN;
if max(demand)<2 || max(assistance)<2, return, end
best = -inf; bestShift = 0;
for shift = -round(.5/dt):round(.5/dt)
    if shift>=0
        a=demand(1:end-shift); b=assistance(1+shift:end);
    else
        a=demand(1-shift:end); b=assistance(1:end+shift);
    end
    a=a-mean(a); b=b-mean(b);
    value=sum(a.*b)/max(norm(a)*norm(b),eps);
    if value>best, best=value; bestShift=shift; end
end
lagMs=bestShift*dt*1000;
end

function index = sustained_onset(signal,threshold)
above = signal>threshold;
index = find(above(1:end-1) & above(2:end),1);
if isempty(index), index=NaN; end
end
