function settings=current_528_settings(subject)
%CURRENT_528_SETTINGS Current-model controller setup and causal input forecast.
% Physical BSE coefficients are read-only. Tuning is fixed across nine cases.
folder=fileparts(mfilename('fullpath'));
settings.controllers=cell(1,2);
% Velocity-form PID: one fixed tuning per activation source, across nine cases.
% RQ1 and RQ2 always use the same closed-loop run and the same tuning.
settings.pid={ ...
    struct('Kp',.83,'Ki',3.4,'Kd',.12,'derivativeTime',.027,'releaseTime',.065), ...
    struct('Kp',.3901,'Ki',7.8257,'Kd',.1786,'derivativeTime',.0197,'releaseTime',.0334)};
% Original movement-PD structure and gains. Only its correction is smoothed.
% P, D, positive correction limit below 1 rad, negative limit, smoothing tau.
settings.humanPD=[5 3 .005 .3 -.001/log(.92)];
settings.actuator=[80 200]; % N m cap and positive rise rate (N m/s).
for source=1:2
    o=struct;
    o.points=40;o.controlPoints=40;o.blocks=8;o.predictionDt=.01;
    if source==2,o.points=120;end % Predict recovery after the 0.4-s control window.
    o.predictionSubsteps=10;o.actualSubsteps=10;o.controlPeriod=.02;
    o.capNm=80;o.rateNmPerSec=200;o.instantRelease=true;
    o.instantActuator=false;o.releaseRateNmPerSec=inf;o.releaseTrackingBandNm=0;
    o.backTorqueMarginNm=.5;o.predictionErrorGuardScale=0;
    o.predictionErrorGuardNm=zeros(o.points,1);o.predictionErrorGuardRelative=[];
    o.predictionErrorGuardFloorNm=0;o.predictionErrorGuardCapNm=0;
    o.angleGovernorUpperRad=0;o.maximumAngleLeadRad=0;
    o.torqueConstraintPriority=true;o.torqueConstraintHorizonSeconds=.4;
    o.angleWeight=20;o.velocityWeight=5;o.backTorqueWeight=300;
    o.assistanceEffortWeight=.1;o.abdomenReferenceWeight=0;
    o.abdomenAllowanceNm=10;o.nonlinearTorqueToleranceNm=0;
    o.inputForecast='arx';o.mapPhaseObserver=false;
    o.activationResidualPredictor=fit_forecast(subject,o.points,o.predictionDt,source);
    predictorFile=fullfile(folder,'data','current_activation_predictor_528.mat');
    if source==1 && isfile(predictorFile)
        data=load(predictorFile);
        modelFile=fullfile(folder,'..','BSE','Direct_Activation_No_Limits','data','subject_528.mat');
        assert(strcmp(data.modelSHA256,nmpc_file_hash(modelFile)), ...
            'The activation predictor belongs to a different BSE input revision.');
        o.activationResidualPredictor=data.predictor;
    end
    o.comparisonVersion='current-dual-MTU-gated-PD-528-v1';
    o.currentCalibratedModel=true;
    settings.controllers{source}=o;
end
settings.protocol=struct('thresholdFraction',.75,'thresholdMuscle','back', ...
    'target','min(unassisted back torque, threshold)','matchedEffort',false, ...
    'appliedTorqueGatesMovementPD',true,'standingPD',true,'angleBoundaries',false);
end

function predictor=fit_forecast(subject,points,predictionDt,source)
% Fit activation residuals only; future experimental inputs are training labels.
% At runtime only past residuals and the trained MAP plan enter the forecast.
lags=0:2:40;X=[];Y=[];
for group=1:numel(subject.sequences)
    q=subject.sequences{group};model=subject.models{source};
    observed=prepare_activation(q,model,source);
    planned=prepare_activation(q,subject.models{2},2);
    residual=observed-planned;residual(q.waiting>.5,:)=0;
    n=numel(q.time);current=(1:2:n)';features=zeros(numel(current),2*numel(lags));
    for k=1:numel(lags)
        history=max(1,current-lags(k));features(:,2*k-1:2*k)=residual(history,:);
    end
    target=zeros(numel(current),2*points);
    for k=1:points
        query=min(n,current+k*predictionDt/q.dt);
        target(:,2*k-1:2*k)=interp1((1:n)',residual,query,'linear');
    end
    X=[X;features];Y=[Y;target]; %#ok<AGROW>
end
scale=max(std(X,1),.002);X=X./scale;
penalty=.01;coefficients=(X'*X+penalty*eye(size(X,2)))\(X'*Y);
predictor=struct('lags',lags,'scale',scale,'coefficients',coefficients, ...
    'points',points,'predictionDt',predictionDt,'penalty',penalty, ...
    'usesOnlyPastActivationResiduals',true,'usesAngleOrTorqueLabels',false, ...
    'calibrationRMSE',sqrt(mean((X*coefficients-Y).^2,'all')), ...
    'baselineMAPResidualRMSE',sqrt(mean(Y.^2,'all')));
end
