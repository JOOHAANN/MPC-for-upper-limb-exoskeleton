% Publish one source-specific NMPC tuning; no PID-effort matching is used.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','NMPC');
old=load(fullfile(demo,'data','comparison_settings_528.mat'));
selected=load(fullfile(folder,'experimental_release28_band20_lead0.2_early1_floor3.5_group1.mat'));
controllers=cell(1,2); controllers{1}=selected.options;
selected=load(fullfile(folder,'map_tracking_group1.mat'));
controllers{2}=selected.options;
controllers{2}.releaseRateNmPerSec=200;
controllers{2}.releaseTrackingBandNm=0;
controllers{2}.maximumAngleLeadRad=.2;
controllers{2}.angleGovernorUpperRad=1.7;
% Pose preservation belongs in the NMPC objective, not a post-QP torque cut.
controllers{1}.angleGovernorUpperRad=0;
controllers{2}.angleGovernorUpperRad=0;
controllers{1}.maximumAngleLeadRad=0;
controllers{2}.maximumAngleLeadRad=0;
controllers{1}.angleWeight=500;
controllers{2}.angleWeight=1000;
controllers{1}.velocityWeight=500;
controllers{2}.velocityWeight=500;
controllers{1}.releaseRateNmPerSec=20;
controllers{1}.releaseTrackingBandNm=40;
controllers{2}.predictionErrorGuardFloorNm=3.5;
controllers{1}.controlPoints=16;
controllers{2}.controlPoints=16;
controllers{2}.points=48;
predictor=controllers{2}.activationResidualPredictor;
originalCoefficients=predictor.coefficients;
predictor.coefficients=zeros(size(originalCoefficients,1),2*controllers{2}.points);
predictor.coefficients(:,1:32)=originalCoefficients;
for point=17:controllers{2}.points
    decay=exp(-(point-16)*controllers{2}.predictionDt/.15);
    predictor.coefficients(:,2*point-1:2*point)=decay*originalCoefficients(:,end-1:end);
end
predictor.points=controllers{2}.points;
controllers{2}.activationResidualPredictor=predictor;
fields={'predictionErrorGuardNm','predictionErrorGuardRelative'};
for k=1:numel(fields)
    values=controllers{2}.(fields{k});
    controllers{2}.(fields{k})=[values(:);repmat(values(end),32,1)];
end
unused={'activationUncertainty','activationEdgeTimesSeconds'};
for source=1:2
    % Resolve ideal motor steps on the same time grid as the actual plant.
    % Prediction durations and the 0.4-s optimized control horizon stay fixed.
    o=controllers{source}; pred=o.activationResidualPredictor;
    oldTime=(0:o.points)*o.predictionDt;
    newPoints=round(o.points*o.predictionDt/.01);
    newTime=(1:newPoints)*.01;
    coefficients=zeros(size(pred.coefficients,1),2*newPoints);
    for channel=1:2
        atZero=zeros(size(pred.coefficients,1),1); atZero(channel)=pred.scale(channel);
        values=[atZero pred.coefficients(:,channel:2:end)];
        coefficients(:,channel:2:end)=interp1(oldTime,values',newTime,'linear')';
    end
    pred.coefficients=coefficients; pred.points=newPoints; pred.predictionDt=.01;
    controllers{source}.activationResidualPredictor=pred;
    controllers{source}.points=newPoints;
    controllers{source}.controlPoints=40;
    controllers{source}.predictionDt=.01;
    controllers{source}.predictionSubsteps=10;
    controllers{source}.instantActuator=false;
    controllers{source}.instantRelease=true;
    controllers{source}.rateNmPerSec=200;
    controllers{source}.releaseRateNmPerSec=inf;
    controllers{source}.releaseTrackingBandNm=0;
    controllers{source}.backTorqueMarginNm=.5;
    controllers{source}.predictionErrorGuardScale=0;
    controllers{source}.predictionErrorGuardFloorNm=0;
    controllers{source}.predictionErrorGuardNm=zeros(controllers{source}.points,1);
    controllers{source}.predictionErrorGuardRelative=[];
    controllers{source}.nonlinearTorqueToleranceNm=0;
    controllers{source}.backTorqueWeight=1800;
    controllers{source}.backTorqueUnderWeight=0;
    controllers{source}.trackingHeadroomNm=0;
    for k=1:numel(unused)
        if isfield(controllers{source},unused{k})
            controllers{source}=rmfield(controllers{source},unused{k});
        end
    end
    controllers{source}.comparisonVersion='clipped-torque-asymmetric-528-v5';
end
selectedPID=load(fullfile(demo,'data','pid_damped_528.mat'));
pid=selectedPID.pid; modelSHA256=old.modelSHA256;
protocol=struct('comparison','shared actuator constraints','effortMatched',false, ...
    'oneNMPCPerSource',true,'conditionsPerSource',9,'sameTracesForRQ1AndRQ2',true, ...
    'thresholdFraction',.75,'pidGainsUnchanged',false,'samePIDAcrossCasesAndRQs',true, ...
    'pidTuning','Damping-focused PID with rise-rate-aware conditional integration', ...
    'actuator','200 N m/s torque increase; instantaneous withdrawal; shared 0-80 N m cap', ...
    'torqueTarget','min(unassisted BSE torque, displayed muscle threshold)', ...
    'tuningScope','all nine recorded calibration scenarios; no independent test set', ...
    'sourceNames',{{'Experimental activation','MAP activation'}});
save(fullfile(demo,'data','threshold_tracking_settings_528.mat'), ...
    'controllers','pid','protocol','modelSHA256','-v7');
rows=struct([]);
for source=1:2
    o=controllers{source};
    row=struct('Source',string(protocol.sourceNames{source}), ...
        'CapNm',o.capNm,'TorqueIncreaseRateNmPerSec',o.rateNmPerSec, ...
        'ControlPeriodSeconds',o.controlPeriod,'PredictionPoints',o.points, ...
        'PredictionDtSeconds',o.predictionDt,'ControlBlocks',o.blocks, ...
        'ControlPredictionPoints',o.controlPoints, ...
        'AngleWeight',o.angleWeight,'VelocityWeight',o.velocityWeight, ...
        'BackTorqueWeight',o.backTorqueWeight,'BackTorqueMarginNm',o.backTorqueMarginNm, ...
        'GuardScale',o.predictionErrorGuardScale,'GuardFloorNm',o.predictionErrorGuardFloorNm, ...
        'UnderTargetWeight',o.backTorqueUnderWeight,'TargetHeadroomNm',o.trackingHeadroomNm, ...
        'AssistanceRegularizationWeight',o.assistanceEffortWeight, ...
        'ReleaseRateNmPerSec',o.releaseRateNmPerSec, ...
        'FastReleaseTrackingBandNm',o.releaseTrackingBandNm, ...
        'AngleGovernorUpperRad',o.angleGovernorUpperRad, ...
        'RecoveryMaximumAngleLeadRad',o.maximumAngleLeadRad, ...
        'RecoveryBackTargetNm',o.recoveryBackTargetNm);
    rows=[rows;row]; %#ok<AGROW>
end
writetable(struct2table(rows),fullfile(demo,'data','threshold_tracking_controller_settings.csv'));
fprintf('Published fixed experimental/MAP controllers with shared 80-N-m equipment limits.\n');
