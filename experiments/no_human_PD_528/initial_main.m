clc
clearvars -except Source_list Group_list Show_figures Reuse_cache Max_samples NMPC_weights Run_label
close all

%% Subject 528: external torque and trajectory control without human posture PD
% This is a separate experiment. The previous comparison is not modified.
% Both controllers act only through exoskeleton torque. The actual BSE and
% every NMPC prediction have zero movement/standing posture-feedback gains.
% PID receives the current task angle; NMPC receives the planned MAP horizon.
% The same PID/NMPC trajectories are used for RQ1, RQ2 and timing evaluation.

%% Select cases and shared equipment
if ~exist('Source_list','var'), Source_list=1:2; end
if ~exist('Group_list','var'), Group_list=1:3; end
if ~exist('Show_figures','var'), Show_figures=true; end
if ~exist('Reuse_cache','var'), Reuse_cache=true; end
if ~exist('Max_samples','var'), Max_samples=inf; end % Development smoke tests only.
if ~exist('NMPC_weights','var'), NMPC_weights=[500 100 20000 100]; end % Angle, torque, exceedance, velocity.
if ~exist('Run_label','var'), Run_label='no_human_PD_528'; end
folder=fileparts(mfilename('fullpath'));
addpath(folder,fullfile(folder,'..','BSE'));
outputFolder=fullfile(folder,'output',Run_label);
if ~isfolder(outputFolder), mkdir(outputFolder); end
subjects=bse_subjects('528'); subject=subjects{1};
calibration=load(fullfile(folder,'data','threshold_tracking_settings_528.mat'));
taskReference=load(fullfile(folder,'data','no_human_PD_task_528.mat'));
assert(strcmp(calibration.modelSHA256,subject.modelSHA256));
assert(strcmp(taskReference.modelSHA256,subject.modelSHA256));
assert(exist('quadprog','file')==2,'Optimization Toolbox is required.');
names={'Experimental activation','MAP activation'};
u_max=80;                        % Shared signed range: -80 to +80 N m.
u_rate_max=200;                  % Rate limit on increasing torque magnitude.
Control_period=.02;              % 50-Hz external PID and NMPC updates.
Threshold_fraction=.75;          % Unassisted, no-human-PD peak for each case.

%% Fixed PID and NMPC settings across all nine cases
pid=calibration.pid;             % Reuse the damped muscle-torque PID gains.
pid.angleKp=100;                 % External actuator feedback, not human PD.
pid.angleKi=8;
pid.angleKd=25;
options=cell(1,2);
for source=1:2
    o=calibration.controllers{source};
    o=prediction_grid(o,20,.025); % 0.5-s forecast with 4 control blocks.
    o.blocks=4; o.capNm=u_max; o.rateNmPerSec=u_rate_max;
    o.controlPeriod=Control_period; o.actualSubsteps=10;
    o.angleWeight=NMPC_weights(1); o.velocityWeight=NMPC_weights(4);
    o.backTorqueWeight=NMPC_weights(2); o.excessWeight=NMPC_weights(3);
    o.effortWeight=.05; o.moveWeight=.2;
    o.backTorqueMarginNm=.5;
    options{source}=o;
end
% One fixed tuning per source is permitted; no load/speed-specific gains.
signature=source_signature(folder);
results=struct('subject',subject.id,'mode','no_human_PD','pid',pid, ...
    'controllers',{options},'modelSHA256',subject.modelSHA256, ...
    'algorithmSHA256',signature,'groups',{cell(1,3)});
results.protocol=struct('humanPostureFeedback',false,'bidirectionalActuator',true, ...
    'lowerTorqueNm',-80,'upperTorqueNm',80,'magnitudeRiseRateNmPerSec',200, ...
    'instantWithdrawal',true,'controlPeriodSeconds',.02,'effortMatched',false, ...
    'thetaReference',taskReference.provenance, ...
    'PID','Reactive torque PID plus current-angle PID through the same external actuator', ...
    'thresholdFraction',Threshold_fraction,'independentValidation',false);

%% Continuous 5/10/15-kg sequences, with one initial state per sequence
for group=Group_list
    sequence=subject.sequences{group};
    if isfinite(Max_samples), sequence=short_sequence(sequence,Max_samples); end
    Ts=sequence.dt;
    results.groups{group}=struct('sequence',sequence,'runs',{cell(1,2)});
    % A frozen, declared MAP task is common to both controllers and sources.
    % It was generated offline from the successful calibrated MAP rollout.
    % The new actual/prediction plant contains no human posture controller.
    planned=bse_inputs(sequence,subject.models{2},2);
    task=taskReference.tasks{group};
    task.theta=task.theta(1:numel(sequence.time));
    task.velocity=task.velocity(1:numel(sequence.time));
    for source=Source_list
        fprintf('\nNo human PD | group %d | %s\n',group,names{source});
        model=subject.models{source};
        plant=bse_without_posture(model,subject.balanceController,sequence.standingAngle);
        inputs=bse_inputs(sequence,model,source); o=options{source};
        baseline=plant.simulate(inputs,sequence.mass,sequence.waiting,sequence.initialState,Ts,10);
        threshold=muscle_threshold(sequence,baseline,Threshold_fraction);
        name=fullfile(outputFolder,sprintf('group%d_source%d.mat',group,source));
        useCache=false;
        if Reuse_cache && isfile(name)
            cache=load(name);
            useCache=strcmp(cache.signature,signature) && isequaln(cache.options,o) ...
                && isequaln(cache.pid,pid) && isequaln(cache.task,task) ...
                && isequaln(cache.sequence,sequence) && isequaln(cache.inputs,inputs) ...
                && strcmp(cache.modelSHA256,subject.modelSHA256);
        end
        if useCache
            pair=cache.pair;
        else
            pidRun=simulate_controller('PID',plant,sequence,inputs,planned,task,threshold,o,pid);
            nmpcRun=simulate_controller('NMPC',plant,sequence,inputs,planned,task,threshold,o,pid);
            pair=struct('baseline',baseline,'pid',pidRun,'nmpc',nmpcRun, ...
                'threshold',threshold,'task',task,'options',o);
            optionsSaved=o; modelSHA256=subject.modelSHA256; %#ok<NASGU>
            cache=struct('signature',signature,'options',optionsSaved,'pid',pid,'task',task, ...
                'sequence',sequence,'inputs',inputs,'pair',pair,'modelSHA256',modelSHA256);
            save(name,'-struct','cache','-v7');
        end
        check_run(pair,Ts);
        results.groups{group}.runs{source}=pair;
    end
end
results.metrics=case_metrics(results,names);
writetable(results.metrics,fullfile(outputFolder,'case_metrics.csv'));
save(fullfile(outputFolder,'comparison_no_human_PD_528.mat'),'results','-v7');
plot_no_posture_comparison(results,outputFolder,Show_figures);
disp(results.metrics(:,{'Source','LoadKg','SpeedCode','Controller', ...
    'ExcessPeakNm','ThetaTaskRMSEDeg','ThetaTaskMaximumErrorDeg','AssistPeakMagnitudeNm'}));
fprintf('\nResults: %s\n',outputFolder);

function run=simulate_controller(kind,plant,q,inputs,planned,task,threshold,o,pid)
n=size(inputs,1); dt=q.dt; stride=round(o.controlPeriod/dt);
[state,first]=plant.initialize(inputs(1,:),q.mass(1),q.waiting(1),q.initialState);
shadow=state;
run.actual=zeros(n,27); run.actual(1,:)=first;
run.command=zeros(n,1); run.solverExitFlag=nan(n,1); run.solveSeconds=zeros(n,1);
run.target=zeros(n,1); run.target(1)=min(first(4),threshold(1));
command=0; warm=zeros(o.blocks,1); memory=zeros(1,5);
solver=optimoptions('quadprog','Algorithm','interior-point-convex','Display','off', ...
    'ConstraintTolerance',1e-8,'OptimalityTolerance',1e-7,'MaxIterations',100);
nextProgress=.25;
for k=2:n
    current=k-1;
    if mod(k-2,stride)==0
        tick=tic;
        if strcmp(kind,'PID')
            [command,memory]=pid_update(run.actual(current,:),state(29), ...
                threshold(current),task.theta(current),task.velocity(current), ...
                q.waiting(current),memory,pid,o);
            flag=1;
        else
            forecast=nmpc_forecast(q,inputs,planned,current,o);
            query=min(n,max(1,current+(1:o.points)'*o.predictionDt/dt));
            reference.theta=interp1((1:n)',task.theta,query,'linear');
            reference.velocity=interp1((1:n)',task.velocity,query,'linear');
            reference.threshold=interp1((1:n)',threshold,query,'previous');
            [command,plan,flag]=nmpc_update(plant,state,shadow,warm,forecast,reference,o,solver);
            grid=(0:o.blocks-1)'*o.points*o.predictionDt/o.blocks;
            warm=interp1(grid,plan,grid+o.controlPeriod,'linear',plan(end));
        end
        run.solveSeconds(k)=toc(tick); run.solverExitFlag(k)=flag;
    end
    applied=motor_command(state(29),command,o.rateNmPerSec*dt,o.capNm);
    [state,run.actual(k,:)]=plant.step(state,inputs(k,:),q.mass(k),q.waiting(k),applied,dt,o.actualSubsteps);
    [shadow,nominal]=plant.step(shadow,inputs(k,:),q.mass(k),q.waiting(k),0,dt,o.actualSubsteps);
    run.command(k)=command; run.target(k)=min(nominal(4),threshold(k));
    if k/n>=nextProgress
        fprintf('  %s: %3.0f%%\n',kind,100*k/n); nextProgress=nextProgress+.25;
    end
end
run.pid=pid; run.inputs=inputs;
end

function [u,m]=pid_update(actual,current,threshold,qref,vref,waiting,m,p,o)
% m = torque integral, angle integral, filtered torque derivative,
% previous measured torque, and threshold-trigger latch.
dt=o.controlPeriod; error=actual(4)-threshold;
d=exp(-dt/p.derivativeTime);
if m(4)==0, m(4)=actual(4); end
m(3)=d*m(3)+(1-d)*(actual(4)-m(4))/dt; m(4)=actual(4);
if waiting, m(1)=0; m(3)=0; m(5)=0; end
if error>0 && ~waiting, m(5)=1; end
qerror=qref-actual(1); verror=vref-actual(2);
trial=m; trial(2)=max(-30,min(30,m(2)+p.angleKi*dt*qerror));
torque=0;
if m(5)
    trial(1)=max(0,min(o.capNm,m(1)+p.Ki*dt*error));
    torque=max(0,p.Kp*error+trial(1)+p.Kd*m(3));
end
raw=torque+p.angleKp*qerror+trial(2)+p.angleKd*verror;
reachable=motor_command(current,raw,o.rateNmPerSec*dt,o.capNm);
% Conditional integration accounts for total-command and rate saturation.
if (raw-reachable)*(trial(1)+trial(2)-m(1)-m(2))<=0
    m(1:2)=trial(1:2);
end
if m(5), torque=max(0,p.Kp*error+m(1)+p.Kd*m(3)); end
u=max(-o.capNm,min(o.capNm,torque+p.angleKp*qerror+m(2)+p.angleKd*verror));
end

function [u,plan,flag]=nmpc_update(plant,state,shadow,warm,forecast,ref,o,solver)
% One Gauss-Newton SQP step, followed by a nonlinear line search.
n=o.blocks; warm=max(-o.capNm,min(o.capNm,warm));
current=state(29); warm(1)=motor_command(current,warm(1),o.rateNmPerSec*o.controlPeriod,o.capNm);
nominal=plant.preview(shadow,zeros(n,1),forecast,o);
base=plant.preview(state,warm,forecast,o);
plan=warm; u=warm(1); flag=-2;
if isempty(nominal) || isempty(base), return, end
ref.torque=min(nominal(:,4),ref.threshold);
residual=tracking_residual(base,warm,current,ref,o);
jac=zeros(numel(residual),n);
for j=1:n
    trial=warm; delta=.5; if trial(j)>o.capNm-.5, delta=-.5; end
    trial(j)=trial(j)+delta;
    predicted=plant.preview(state,trial,forecast,o);
    if isempty(predicted), continue, end
    jac(:,j)=(tracking_residual(predicted,trial,current,ref,o)-residual)/delta;
end
H=jac'*jac+1e-6*eye(n); f=jac'*residual-H*warm;
lower=-o.capNm*ones(n,1); upper=o.capNm*ones(n,1);
% Requests can withdraw immediately. The plant enforces the magnitude-rise
% limit at every integration interval, including zero-crossing reversals.
lower(1)=-max(-current,0)-o.rateNmPerSec*o.controlPeriod;
upper(1)=max(current,0)+o.rateNmPerSec*o.controlPeriod;
lower(1)=max(-o.capNm,lower(1)); upper(1)=min(o.capNm,upper(1));
[candidate,~,exitFlag]=quadprog(H,f,[],[],[],[],lower,upper,[],solver);
if exitFlag<=0 || isempty(candidate) || any(~isfinite(candidate)), return, end
best=residual'*residual;
for fraction=[1 .5 .25 .125]
    proposed=warm+fraction*(candidate-warm);
    predicted=plant.preview(state,proposed,forecast,o);
    if isempty(predicted), continue, end
    r=tracking_residual(predicted,proposed,current,ref,o);
    if r'*r<best
        best=r'*r; plan=proposed;
    end
end
u=max(lower(1),min(upper(1),plan(1))); plan(1)=u; flag=1;
end

function r=tracking_residual(predicted,controls,current,ref,o)
% Torque ceiling is soft because a single finite actuator may not satisfy
% both targets simultaneously. No actual force, angle or plotted output is clipped.
r=[sqrt(o.angleWeight)*(predicted(:,1)-ref.theta)/.1; ...
    sqrt(o.velocityWeight)*(predicted(:,2)-ref.velocity); ...
    sqrt(o.backTorqueWeight)*(predicted(:,4)-ref.torque)/30; ...
    sqrt(o.excessWeight)*max(predicted(:,4)-ref.threshold+o.backTorqueMarginNm,0)/30; ...
    sqrt(o.effortWeight)*controls/80; ...
    sqrt(o.moveWeight)*diff([current;controls])/20];
end

function u=motor_command(previous,requested,increment,cap)
requested=max(-cap,min(cap,requested));
if requested*previous<0, previous=0; end
u=sign(requested)*min(abs(requested),abs(previous)+increment);
end

function threshold=muscle_threshold(q,baseline,fraction)
threshold=nan(size(q.time)); limits=[];
for k=1:numel(q.segments)
    seg=q.segments{k}; ix=seg.startIndex:min(seg.endIndex,numel(q.time));
    if isempty(ix), continue, end
    limit=fraction*max(baseline(ix,4)); threshold(ix)=limit; limits(end+1)=limit; %#ok<AGROW>
end
if isempty(limits), limits=fraction*max(baseline(:,4)); end
threshold(isnan(threshold))=min(limits);
end

function check_run(pair,dt)
for data={pair.baseline,pair.pid.actual,pair.nmpc.actual}
    a=data{1}; assert(all(isfinite(a),'all'));
    assert(max(abs(a(:,22:23)),[],'all')<1e-10,'Human posture feedback remains active.');
    assert(max(abs(a(:,3)-a(:,4)+a(:,5)))<1e-6);
end
for data={pair.pid,pair.nmpc}
    r=data{1}; u=r.actual(:,25);
    assert(max(abs(u))<=80+1e-8);
    for k=2:numel(u)
        assert(abs(u(k)-motor_command(u(k-1),r.command(k),200*dt,80))<1e-8);
    end
end
end

function t=case_metrics(results,names)
rows=struct([]);
for group=1:3
    if isempty(results.groups{group}), continue, end
    q=results.groups{group}.sequence;
    for source=1:2
        p=results.groups{group}.runs{source}; if isempty(p), continue, end
        actual={p.baseline,p.pid.actual,p.nmpc.actual}; labels={'Unassisted','PID','NMPC'};
        for stage=1:numel(q.segments)
            seg=q.segments{stage}; ix=seg.startIndex:min(seg.endIndex,numel(q.time));
            if isempty(ix), continue, end
            for method=1:3
                a=actual{method}; excess=max(a(ix,4)-p.threshold(ix),0);
                error=a(ix,1)-p.task.theta(ix);
                context=max(1,ix(1)-40):min(numel(q.time),ix(end)+80);
                demand=max(p.baseline(context,4)-p.threshold(context),0);
                [lag,corr]=timing_lag(demand,max(a(context,25),0),q.dt);
                row=struct('Source',string(names{source}),'Group',group,'LoadKg',seg.loadKg, ...
                    'SpeedCode',seg.speedCode,'Controller',string(labels{method}), ...
                    'ThresholdNm',p.threshold(ix(1)),'ExcessPeakNm',max(excess), ...
                    'ExcessAreaNms',sum(excess)*q.dt,'ExcessDurationSeconds',sum(excess>.1)*q.dt, ...
                    'ThetaTaskRMSEDeg',rad2deg(sqrt(mean(error.^2))), ...
                    'ThetaTaskMaximumErrorDeg',rad2deg(max(abs(error))), ...
                    'ThetaFromUnassistedRMSEDeg',rad2deg(sqrt(mean((a(ix,1)-p.baseline(ix,1)).^2))), ...
                    'MinimumThetaRad',min(a(context,1)),'MaximumThetaRad',max(a(context,1)), ...
                    'AssistPeakMagnitudeNm',max(abs(a(context,25))), ...
                    'AssistanceMinimumNm',min(a(context,25)),'AssistanceMaximumNm',max(a(context,25)), ...
                    'SquaredTorqueEffortNm2s',sum(a(context,25).^2)*q.dt, ...
                    'PositiveDemandLagMs',lag,'PositiveDemandCorrelation',corr, ...
                    'AbdomenPeakNm',max(a(context,5)),'SoftBoundaryPeakMagnitudeNm',max(abs(a(context,26))));
                rows=[rows;row]; %#ok<AGROW>
            end
        end
    end
end
t=struct2table(rows);
end

function [lag,best]=timing_lag(demand,assistance,dt)
lag=nan; best=nan;
if max(demand)<2 || max(assistance)<2, return, end
best=-inf;
for shift=-round(.5/dt):round(.5/dt)
    if shift>=0, a=demand(1:end-shift); b=assistance(1+shift:end);
    else, a=demand(1-shift:end); b=assistance(1:end+shift); end
    a=a-mean(a); b=b-mean(b); value=sum(a.*b)/max(norm(a)*norm(b),eps);
    if value>best, best=value; lag=shift*dt*1000; end
end
end

function o=prediction_grid(o,points,dt)
predictor=o.activationResidualPredictor;
old=(0:predictor.points)*predictor.predictionDt;
coefficients=zeros(size(predictor.coefficients,1),2*points);
for channel=1:2
    atZero=zeros(size(coefficients,1),1); atZero(channel)=predictor.scale(channel);
    values=[atZero predictor.coefficients(:,channel:2:end)];
    coefficients(:,channel:2:end)=interp1(old,values',min(old(end),(1:points)*dt),'linear')';
end
predictor.coefficients=coefficients; predictor.points=points; predictor.predictionDt=dt;
o.activationResidualPredictor=predictor; o.points=points; o.predictionDt=dt;
o.predictionSubsteps=round(dt/.001);
end

function q=short_sequence(q,count)
count=min(count,numel(q.time)); fields=fieldnames(q);
originalCount=numel(q.time);
for k=1:numel(fields)
    value=q.(fields{k});
    if isnumeric(value) || islogical(value)
        if size(value,1)==originalCount, q.(fields{k})=value(1:count,:); end
    end
end
end

function value=source_signature(folder)
files={mfilename('fullpath'),fullfile(folder,'bse_without_posture.m'), ...
    fullfile(folder,'private','native','bse_no_posture_mex.c'),fullfile(folder,'nmpc_forecast.m')};
engine=java.security.MessageDigest.getInstance('SHA-256');
for k=1:numel(files)
    name=files{k}; if ~isfile(name), name=[name '.m']; end
    f=fopen(name,'rb'); bytes=fread(f,inf,'*uint8'); fclose(f); engine.update(bytes);
end
value=lower(reshape(dec2hex(typecast(engine.digest(),'uint8'),2)',1,[]));
end
