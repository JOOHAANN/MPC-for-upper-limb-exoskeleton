clc
clearvars -except Subject_list Max_evaluations

%% Refit the direct-activation model, without movement PD or angle limits
% The trained parameter files are the starting point for a repeat fit.
% All 24 physical coefficients can change. Data-derived MAP scaling is fixed.
% Each source keeps one fixed vector for all nine cases.
% Existing data files are backed up before replacement.
if ~exist('Subject_list','var'),Subject_list={'528'};end
if ~exist('Max_evaluations','var'),Max_evaluations=1200;end
folder=fileparts(mfilename('fullpath'));addpath(folder);
output=fullfile(folder,'output','retraining');if ~isfolder(output),mkdir(output);end
assert(exist('lsqnonlin','file')==2,'Optimization Toolbox is required for training.');
options=optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',Max_evaluations, ...
    'FiniteDifferenceStepSize',1e-4,'FunctionTolerance',1e-4,'StepTolerance',1e-4);
for subjectIndex=1:numel(Subject_list)
    file=fullfile(folder,'data',['subject_' Subject_list{subjectIndex} '.mat']);
    copyfile(file,fullfile(output,['subject_' Subject_list{subjectIndex} '_before.mat']));
    data=load(file);subject=data.subject;
    for source=1:2
        model=subject.models{source};assert(numel(model.parameters)==24);
        x0=parameter_scale(model.parameters,'encode');
        objective=@(x)calibration_error(parameter_scale(x,'decode'),subject,source);
        [x,cost,~,exitFlag,solver]=lsqnonlin(objective,x0,zeros(1,24),ones(1,24),options);
        fitted=parameter_scale(x,'decode');
        subject.models{source}.parameters=fitted(1:24);
        subject.models{source}.abdomenScale=1;
        report=struct('source',source,'cost',cost,'exitFlag',exitFlag,'solver',solver, ...
            'parameters',subject.models{source}.parameters,'abdomenScale',1);
        save(fullfile(output,sprintf('%s_source%d_fit.mat',subject.shortName,source)),'report');
    end
    subject.parameterPolicy='New subject/source-specific calibration over nine recorded bouts; no movement feedback or angle boundaries.';
    subject.calibrationProtocol=struct('movementPD',false,'standingPD',true, ...
        'angleBoundaries',false,'assistanceNm',0,'numberOfPhysicalCoefficients',24, ...
        'numberOfInputGains',0,'numberOfRecordedBouts',9,'outerStepSeconds',.01, ...
        'internalStepSeconds',.001,'independentValidation',false, ...
        'referenceUse','Training loss and evaluation only; never runtime movement feedback.', ...
        'status','calibrated');
    save(file,'subject','-v7');
end
fprintf('Training finished. Run run_dual_MTU_BSE.m to regenerate all figures and metrics.\n');

function residual=calibration_error(fitted,subject,source)
% Reference values appear only here, not in dual_MTU_model.
residual=[];thetaCorrelations=[];torqueCorrelations=[];
p=fitted(1:24);model=subject.models{source};
model.parameters=p;model.abdomenScale=1;
w=calibration_settings(model,source);
for group=1:numel(subject.sequences)
    q=subject.sequences{group};inputs=prepare_activation(q,model,source);
    plant=dual_MTU_model(p,subject.balanceController,q.standingAngle,true);
    out=plant.simulate(inputs,q.mass,q.waiting,q.initialState,q.dt,10);
    for stage=1:numel(q.segments)
        seg=q.segments{stage};available=seg.startIndex:seg.endIndex;
        ix=available(unique(round(linspace(1,numel(available),80))));
        predicted=out(ix,1);predicted(~isfinite(predicted))=10;
        reference=q.referenceTheta(ix);scale=sqrt(numel(ix));
        residual=[residual;w.thetaPoint*(predicted-reference)/(.3*scale); ...
            w.thetaShape*(standardize(predicted,.12)-standardize(reference,.12))/scale]; %#ok<AGROW>
        thetaCorrelations(end+1)=safe_correlation(predicted,reference); %#ok<AGROW>
        it=ix(q.torqueValid(ix) & isfinite(q.referenceTorque(ix)));
        if ~isempty(it)
            predictedTorque=out(it,3);predictedTorque(~isfinite(predictedTorque))=3000;
            referenceTorque=q.referenceTorque(it);scale=sqrt(numel(it));
            residual=[residual;w.torquePoint*(predictedTorque-referenceTorque)/(50*scale); ...
                w.torqueShape*(standardize(predictedTorque,15)-standardize(referenceTorque,15))/scale]; %#ok<AGROW>
            torqueCorrelations(end+1)=safe_correlation(predictedTorque,referenceTorque); %#ok<AGROW>
        end
        back=out(ix,4);abdomen=out(ix,5);
        back(~isfinite(back))=1000;abdomen(~isfinite(abdomen))=1000;
        ratio=mean(abdomen)/max(mean(back),1);
        residual=[residual;w.abdomenRatioWeight*max(ratio-w.abdomenMeanRatioSoftPreference,0); ...
            w.peakRatioWeight*max(max(abdomen)/max(max(back),1)-.8,0)]; %#ok<AGROW>
        if w.extremumWeight>0
            middle=floor(numel(ix)/2);pred=out(ix,1);pred(~isfinite(pred))=10;ref=q.referenceTheta(ix);
            residual=[residual;w.extremumWeight/w.extremumScaleRad*[ ...
                min(pred(1:middle))-min(ref(1:middle)); ...
                min(pred(middle+1:end))-min(ref(middle+1:end));max(max(pred)-max(ref),0)]]; %#ok<AGROW>
        end
        if w.standingRecoverySeconds>0
            post=ix(end)+1:min(numel(q.time),ix(end)+round(w.standingRecoverySeconds/q.dt));
            post=post(q.waiting(post)>.5);
            if ~isempty(post)
                angle=out(post,1);velocity=out(post,2);
                angle(~isfinite(angle))=10;velocity(~isfinite(velocity))=10;
                residual=[residual;(angle-q.standingAngle)/(w.recoveryAngleScaleRad*sqrt(numel(post))); ...
                    .2*velocity/sqrt(numel(post))]; %#ok<AGROW>
            end
        end
    end
    if isfinite(w.restingAbdomenPreferenceNm)
        resting=out(1,5);if ~isfinite(resting),resting=1000;end
        residual=[residual;max(resting-w.restingAbdomenPreferenceNm,0)/ ...
            (w.restingAbdomenScaleNm*sqrt(numel(subject.sequences)))]; %#ok<AGROW>
    end
end
% Broad engineering priors, not independently measured anatomical values.
residual=[residual;.08*log((p(1)/p(17))/6000);.08*log((p(9)/p(18))/1500)];
residual=[residual;w.thetaGoalWeight*max(w.meanThetaTarget-mean(thetaCorrelations),0); ...
    w.torqueGoalWeight*max(w.meanTorqueTarget-mean(torqueCorrelations),0)];
end

function z=standardize(a,floorValue)
z=(a-mean(a))/max(std(a,1),floorValue);
end
function value=safe_correlation(a,b)
value=0;if std(a)<1e-10 || std(b)<1e-10,return,end
c=corrcoef(a,b);value=c(1,2);
end
function w=calibration_settings(model,source)
w=struct('thetaPoint',1.5,'thetaShape',2.5,'torquePoint',1.2,'torqueShape',2.5, ...
    'abdomenRatioWeight',8,'abdomenMeanRatioSoftPreference',.25,'peakRatioWeight',1, ...
    'meanThetaTarget',.73,'meanTorqueTarget',.83,'thetaGoalWeight',70,'torqueGoalWeight',70, ...
    'extremumWeight',.5,'extremumScaleRad',.2,'restingAbdomenPreferenceNm',5, ...
    'restingAbdomenScaleNm',5,'standingRecoverySeconds',1.2,'recoveryAngleScaleRad',.2);
if source==2,w.meanThetaTarget=.66;w.meanTorqueTarget=.79;end
if isfield(model,'fittingWeights')
    fields=fieldnames(model.fittingWeights);
    for k=1:numel(fields),w.(fields{k})=model.fittingWeights.(fields{k});end
end
end

function value=parameter_scale(value,direction)
% Optimize slack tendon lengths; decode them to reference MT lengths.
low=[60 .15 .45 .5 2 .4 .02 .005 5 .4 .25 .45 .345 0 ...
    15 .20 .020 .035 .025 .025 2 .02 .005 .10];
high=[1600 1.1 1.45 60 15 3 .10 .20 650 2 .95 1.45 .85 1.1 ...
    45 .35 .060 .10 .35 .40 15 .10 .20 .35];
value=value(:)';
if strcmp(direction,'encode')
    w=sin(value(24));
    value(19)=value(19)-value(17)/value(2)*sqrt(value(3)^2-w^2);
    value(20)=value(20)-value(18)/value(11)*sqrt(value(12)^2-w^2);
    value=max(1e-6,min(1-1e-6,(value-low)./(high-low)));
else
    value=low+value.*(high-low);w=sin(value(24));
    value(19)=value(19)+value(17)/value(2)*sqrt(value(3)^2-w^2);
    value(20)=value(20)+value(18)/value(11)*sqrt(value(12)^2-w^2);
end
end
