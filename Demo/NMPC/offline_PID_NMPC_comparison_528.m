clc
clearvars -except Source_list Show_figures Reuse_cache
close all

%% PID versus NMPC using the latest subject-528 dual-MTU BSE
% Run this script directly. Both activation sources use all nine conditions.
% No movement PD at zero assistance. Standing/waiting feedback is retained.
% PID and NMPC share the actuator cap and positive torque-rate limit.
% RQ1 and RQ2 use the same saved closed-loop trajectories.

%% User settings
if ~exist('Source_list','var'), Source_list=1:2; end   % 1: experimental; 2: MAP
if ~exist('Show_figures','var'), Show_figures=true; end
if ~exist('Reuse_cache','var'), Reuse_cache=true; end
assert(~isempty(Source_list) && all(ismember(Source_list,1:2)), ...
    'Source_list must contain 1 (Experimental), 2 (MAP), or both.');
sources=unique(Source_list(:)');
showFigures=Show_figures;
reuseCache=Reuse_cache;

%% Load the current calibrated model and fixed controller settings
folder=fileparts(mfilename('fullpath'));addpath(folder,fileparts(folder));
bseFolder=fullfile(folder,'..','BSE','Direct_Activation_No_Limits');addpath(bseFolder);
file=fullfile(bseFolder,'data','subject_528.mat');data=load(file);subject=data.subject;
assert(strcmp(subject.calibrationProtocol.status,'calibrated'));
assert(isfield(subject,'MAPScaling'),'Use the latest legacy-input/load-scaled BSE package.');
assert(exist('quadprog','file')==2,'Optimization Toolbox is required.');
settings=current_528_settings(subject);
label='current_BSE_528';
output=fullfile(folder,'output',label);if ~isfolder(output),mkdir(output);end
signature=struct('model',nmpc_file_hash(file), ...
    'plant',nmpc_file_hash(fullfile(folder,'private','native','current_mtu_control.c')), ...
    'mechanics',nmpc_file_hash(fullfile(bseFolder,'private','native','mtu_mechanics.h')), ...
    'nmpc',nmpc_file_hash(fullfile(folder,'simulate_threshold_nmpc.m')), ...
    'forecast',nmpc_file_hash(fullfile(folder,'nmpc_forecast.m')));
pidAlgorithm=nmpc_file_hash(fullfile(folder,'simulate_threshold_pid.m'));
names={'Experimental activation','MAP activation'};
results=struct('subject',subject.id,'mode','threshold_tracking','pid',{settings.pid}, ...
    'protocol',settings.protocol,'controllers',{settings.controllers},'humanPD',settings.humanPD, ...
    'modelSHA256',signature.model,'signature',signature,'groups',{cell(1,3)});
rows=struct([]);
results.models=subject.models;results.balanceController=subject.balanceController;
results.MAPScaling=subject.MAPScaling;
%% Replay the three continuous lifting sequences
for group=1:3
    q=subject.sequences{group};results.groups{group}=struct('sequence',q,'runs',{cell(1,2)});
    for source=sources
        model=subject.models{source};options=settings.controllers{source};pidSettings=settings.pid{source};
        plant=current_bse_plant(model,subject.balanceController,q.standingAngle,settings.humanPD,settings.actuator);
        inputs=prepare_activation(q,model,source);planned=prepare_activation(q,subject.models{2},2);
        baseline=plant.simulate(inputs,q.mass,q.waiting,q.initialState,q.dt,options.actualSubsteps);
        direct=dual_MTU_model(model.parameters,subject.balanceController,q.standingAngle,true);
        expected=direct.simulate(inputs,q.mass,q.waiting,q.initialState,q.dt,options.actualSubsteps);
        assert(max(abs(baseline-expected),[],'all')<1e-7, ...
            'Current controller plant no longer matches the calibrated unassisted BSE.');
        threshold=muscle_threshold(q,baseline);
        cacheFile=fullfile(output,sprintf('group%d_source%d.mat',group,source));cached=false;
        if reuseCache && isfile(cacheFile)
            cache=load(cacheFile);
            cached=isequaln(cache.signature,signature) && isequaln(cache.options,options) ...
                && isequaln(cache.humanPD,settings.humanPD) ...
                && isequaln(cache.inputs,inputs) && isequaln(cache.planned,planned);
        end
        fprintf('\nCurrent BSE + gated PD | group %d | %s | cache %d\n',group,names{source},cached);
        if cached
            nmpcRun=cache.pair.nmpc;
        else
            nmpcRun=simulate_threshold_nmpc(plant,q,inputs,planned,threshold,options);
            % Preserve the expensive NMPC replay independently of PID tuning.
            pair=struct('baseline',baseline,'nmpc',nmpcRun,'pid',[], ...
                'threshold',threshold,'options',options);
            humanPD=settings.humanPD;pid=[]; %#ok<NASGU>
            save(cacheFile,'pair','signature','options','humanPD','pid','pidAlgorithm','inputs','planned','-v7');
        end
        if cached && isequaln(cache.pid,pidSettings) && isfield(cache,'pidAlgorithm') && strcmp(cache.pidAlgorithm,pidAlgorithm)
            pidRun=cache.pair.pid;
        else
            pidRun=simulate_threshold_pid(plant,q,inputs,threshold,options,pidSettings);
        end
        verify_run(pidRun,q,inputs,settings.actuator);verify_run(nmpcRun,q,inputs,settings.actuator);
        assert(max(abs(nmpcRun.nominal(:,1)-baseline(:,1)))<1e-8);
        pair=struct('baseline',baseline,'pid',pidRun,'nmpc',nmpcRun,'threshold',threshold,'options',options);
        humanPD=settings.humanPD;pid=pidSettings; %#ok<NASGU>
        save(cacheFile,'pair','signature','options','humanPD','pid','pidAlgorithm','inputs','planned','-v7');
        results.groups{group}.runs{source}=pair;
        rows=[rows;pid_nmpc_metrics(q,pair,group,source,names{source})]; %#ok<AGROW>
    end
end
results.metrics=struct2table(rows);results.metrics=rate_and_timing_metrics(results);
results.summary=summarize(results.metrics,sources,names);
writetable(results.metrics,fullfile(output,'comparison_metrics_528.csv'));
writetable(results.summary,fullfile(output,'comparison_summary_528.csv'));
save(fullfile(output,'current_BSE_NMPC_528.mat'),'results','settings','-v7');
%% Export the original RQ1/RQ2 plots and current-model checks
if isequal(sources,1:2)
    save(fullfile(folder,'data','current_controller_528.mat'),'settings','signature','-v7');
end
plot_pid_nmpc_comparison(results,output,showFigures);
if isequal(sources,1:2)
    subject.nmpc.capNm=settings.actuator(1);
    runs=cellfun(@(p)p.nmpc,results.groups{1}.runs,'UniformOutput',false);
    plot_nmpc_results(subject,subject.sequences{1},runs,output,showFigures);
    plot_current_PD(results,output,showFigures);
    [~,settlingSummary]=pid_settling_report(results,output,showFigures);
    write_controller_settings(settings,fullfile(folder,'data','current_controller_settings_528.csv'));
    write_current_audit(results,settings,settlingSummary,fullfile(output,'release_audit.json'));
end
disp(results.summary);
fprintf('\nFigures and metrics: %s\n',output);

%% Local checks and summary calculations
function threshold=muscle_threshold(q,baseline)
threshold=zeros(size(q.time));limits=zeros(1,3);
for k=1:3
    seg=q.segments{k};ix=seg.startIndex:seg.endIndex;
    limits(k)=.75*max(baseline(ix,4));threshold(ix)=limits(k);
end
threshold(q.waiting>.5)=min(limits);
end

function verify_run(run,q,inputs,actuator)
a=run.actual;u=a(:,25);moving=q.waiting<.5;off=moving & abs(u)<=1e-9;
assert(all(isfinite(a),'all') && min(u)>=-1e-9 && max(u)<=actuator(1)+1e-8);
assert(max(diff(u))/q.dt<=actuator(2)+1e-6);
assert(all(a(:,26)==0) && all(a(moving,20:22)==0,'all'));
assert(max(abs(a(off,18:19)-inputs(off,:)),[],'all')<1e-12, ...
    'Movement feedback remained active without assistance.');
assert(all(a(off,24)==0));
assert(max(abs(a(:,3)-a(:,4)+a(:,5)))<1e-6);
assert(max(abs(a(:,16)))<1e-4 && max(abs(a(:,17)))<1e-5);
end

function summary=summarize(metrics,sources,names)
rows=struct([]);
for source=sources
    p=metrics(metrics.Source==string(names{source}) & metrics.Controller=="PID",:);
    n=metrics(metrics.Source==string(names{source}) & metrics.Controller=="NMPC",:);
    row=struct('Source',string(names{source}),'Conditions',height(p), ...
        'PIDMaximumExcessNm',max(p.ExcessPeakNm),'NMPCMaximumExcessNm',max(n.ExcessPeakNm), ...
        'PIDTotalExcessAreaNms',sum(p.ExcessAreaNms),'NMPCTotalExcessAreaNms',sum(n.ExcessAreaNms), ...
        'PIDMeanThetaRMSEDeg',mean(p.ThetaDeviationRMSEDeg),'NMPCMeanThetaRMSEDeg',mean(n.ThetaDeviationRMSEDeg), ...
        'PIDMeanDemandLagMs',mean(p.DemandLagMs,'omitnan'),'NMPCMeanDemandLagMs',mean(n.DemandLagMs,'omitnan'));
    rows=[rows;row]; %#ok<AGROW>
end
summary=struct2table(rows);
end

function metrics=rate_and_timing_metrics(results)
metrics=results.metrics;n=height(metrics);
metrics.MaximumAssistRiseRateNmPerSec=zeros(n,1);
metrics.MaximumAssistWithdrawalRateNmPerSec=zeros(n,1);
metrics.SolveMeanMs=zeros(n,1);metrics.SolveP95Ms=zeros(n,1);metrics.SolveMaxMs=zeros(n,1);
metrics.UpdatesAbove20ms=zeros(n,1);metrics.TotalUpdates=zeros(n,1);
for row=1:n
    group=metrics.Group(row);source=metrics.SourceIndex(row);stage=metrics.Stage(row);
    q=results.groups{group}.sequence;pair=results.groups{group}.runs{source};
    cuts=[1 round((q.segments{1}.endIndex+q.segments{2}.startIndex)/2)+1 ...
        round((q.segments{2}.endIndex+q.segments{3}.startIndex)/2)+1 numel(q.time)+1];
    ix=cuts(stage):cuts(stage+1)-1;
    if metrics.Controller(row)=="Unassisted",continue,end
    method=lower(char(metrics.Controller(row)));run=pair.(method);
    rate=diff(run.actual(ix,25))/q.dt;
    metrics.MaximumAssistRiseRateNmPerSec(row)=max([0;rate]);
    metrics.MaximumAssistWithdrawalRateNmPerSec(row)=max([0;-rate]);
    seconds=run.solveSeconds(ix);seconds=seconds(seconds>0);
    if isempty(seconds),continue,end
    metrics.SolveMeanMs(row)=1000*mean(seconds);metrics.SolveP95Ms(row)=1000*prctile(seconds,95);
    metrics.SolveMaxMs(row)=1000*max(seconds);metrics.UpdatesAbove20ms(row)=sum(seconds>.02);
    metrics.TotalUpdates(row)=numel(seconds);
end
end

function write_controller_settings(settings,file)
% Export what is actually used for each activation source.
names={'Experimental activation','MAP activation'};rows=struct([]);
for source=1:2
    o=settings.controllers{source};p=settings.pid{source};
    row=struct('Source',string(names{source}),'OuterSampleIntervalSeconds',.01, ...
        'InternalStepSeconds',.001,'ControllerUpdateSeconds',o.controlPeriod, ...
        'PredictionPoints',o.points,'PredictionStepSeconds',o.predictionDt, ...
        'PredictionHorizonSeconds',o.points*o.predictionDt,'ControlledPoints',o.controlPoints, ...
        'ControlBlocks',o.blocks,'ControlWindowSeconds',o.controlPoints*o.predictionDt, ...
        'ZeroCommandRecoveryTailSeconds',(o.points-o.controlPoints)*o.predictionDt, ...
        'AssistanceMinimumNm',0,'AssistanceMaximumNm',o.capNm, ...
        'PositiveRiseLimitNmPerSec',o.rateNmPerSec,'WithdrawalRateLimitNmPerSec',inf, ...
        'ThresholdFraction',settings.protocol.thresholdFraction,'FeasibilityReserveNm',o.backTorqueMarginNm, ...
        'MovementPD_P',settings.humanPD(1),'MovementPD_D',settings.humanPD(2), ...
        'PDOnlyWithAppliedAssistance',true,'PID_Kp',p.Kp,'PID_Ki',p.Ki,'PID_Kd',p.Kd, ...
        'PID_DerivativeFilterSeconds',p.derivativeTime,'PID_WithdrawalFilterSeconds',p.releaseTime);
    rows=[rows;row]; %#ok<AGROW>
end
writetable(struct2table(rows),file);
end

function write_current_audit(results,settings,settlingSummary,file)
% All values below are recalculated from the current saved simulations.
audit=struct('subject',results.subject,'modelSHA256',results.modelSHA256, ...
    'scope','Offline calibration-set replay; not a hardware or unseen-condition validation.', ...
    'controller','Reactive velocity-form PID; one fixed tuning per activation source.', ...
    'thresholdFraction',settings.protocol.thresholdFraction,'actuatorLimits',settings.actuator, ...
    'backActivationTimeSeconds',[results.models{1}.parameters(8),results.models{2}.parameters(8)], ...
    'pid',{settings.pid}, ...
    'bandDefinition','Post-peak entry held >=0.15 s within +/-max(5 N m, 5% threshold), for unassisted demand windows >=0.30 s. Stayed means staying inside after first entry. Settled means staying inside after the final re-entry until demand ends. This is not an LTI step-response settling time.');
names={'Experimental activation','MAP activation'};sourceRows=cell(1,2);
for source=1:2
    m=results.metrics(results.metrics.SourceIndex==source,:);
    p=m(m.Controller=="PID",:);n=m(m.Controller=="NMPC",:);
    b=m(m.Controller=="Unassisted",:);
    row=struct('source',names{source},'conditions',height(p), ...
        'pidMaximumExcessNm',max(p.ExcessPeakNm),'nmpcMaximumExcessNm',max(n.ExcessPeakNm), ...
        'pidTotalExcessAreaNms',sum(p.ExcessAreaNms),'nmpcTotalExcessAreaNms',sum(n.ExcessAreaNms), ...
        'unassistedTotalExcessAreaNms',sum(b.ExcessAreaNms), ...
        'pidMeanThetaDeviationRMSEDeg',mean(p.ThetaDeviationRMSEDeg), ...
        'nmpcMeanThetaDeviationRMSEDeg',mean(n.ThetaDeviationRMSEDeg), ...
        'maximumPIDAssistanceNm',max(p.AssistPeakNm),'maximumNMPCAssistanceNm',max(n.AssistPeakNm), ...
        'nmpcUpdatesAbove20ms',sum(n.UpdatesAbove20ms),'nmpcTotalUpdates',sum(n.TotalUpdates));
    row.bandEntry=table2struct(settlingSummary(settlingSummary.Source==string(names{source}),:));
    sourceRows{source}=row;
end
audit.sources=sourceRows;
fid=fopen(file,'w');assert(fid>=0,'Could not write current audit.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(audit,'PrettyPrint',true));
end
