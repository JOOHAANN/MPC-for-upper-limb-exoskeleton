function test_pid_nmpc_comparison
%TEST_PID_NMPC_COMPARISON One tuning per source, shared hardware, nine cases.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'Demo','NMPC');
addpath(folder,fullfile(root,'Demo','BSE'));
data=load(fullfile(folder,'output','pid_nmpc_528','threshold_tracking', ...
    'pid_nmpc_comparison_528.mat'));
r=data.results;
assert(strcmp(r.mode,'threshold_tracking'));
assert(height(r.metrics)==54,'Expected 9 cases x 2 sources x 3 controllers.');
settings=load(fullfile(folder,'data','threshold_tracking_settings_528.mat'));
assert(isequaln(settings.pid,r.pid),'The selected PID must be identical in all comparisons.');
assert(settings.protocol.samePIDAcrossCasesAndRQs);
subjects=bse_subjects('528'); subject=subjects{1};
assert(strcmp(r.modelSHA256,subject.modelSHA256));
poseMaximum=zeros(3,2);
torqueExcess=zeros(3,2);
for source=1:2
    for group=1:3
        sequence=r.groups{group}.sequence;
        pair=r.groups{group}.runs{source};
        assert(isequaln(pair.options,settings.controllers{source}));
        assert(isequaln(pair.pid.pid,r.pid));
        assert(pair.options.capNm==80 && pair.options.rateNmPerSec==200);
        assert(pair.options.instantRelease && ~pair.options.instantActuator);
        assert(pair.options.backTorqueMarginNm==.5);
        assert(pair.options.predictionErrorGuardFloorNm==0 && pair.options.predictionErrorGuardScale==0);
        assert(pair.options.controlPeriod==.02);
        assert(isinf(pair.options.releaseRateNmPerSec));
        assert(max(abs(pair.nmpc.nominal-pair.baseline),[],'all')<1e-7);
        target=min(pair.baseline(:,4),pair.threshold);
        assert(max(abs(pair.nmpc.target-target))<1e-8,'Clipped reference was shifted or replaced.');
        below=pair.baseline(:,4)<pair.threshold;
        assert(max(abs(pair.nmpc.target(below)-pair.baseline(below,4)))<1e-8);
        torqueExcess(group,source)=max(pair.nmpc.actual(:,4)-pair.threshold);
        poseMaximum(group,source)=max(pair.nmpc.actual(:,1));
        for actual={pair.pid.actual,pair.nmpc.actual}
            a=actual{1}; u=a(:,25);
            assert(all(isfinite(a),'all'));
            assert(min(u)>=-1e-10 && max(u)<=80+1e-8);
            assert(max(abs(a(:,3)-a(:,4)+a(:,5)))<1e-6);
        end
        for controlled={pair.pid,pair.nmpc}
            c=controlled{1}; u=c.actual(:,25);
            expected=min(c.command(2:end),u(1:end-1)+200*sequence.dt);
            assert(max(abs(u(2:end)-expected))<1e-8);
            assert(max(diff(u))/sequence.dt<=200+1e-6);
        end
        for stage=1:3
            segment=sequence.segments{stage}; ix=segment.startIndex:segment.endIndex;
            threshold=.75*max(pair.baseline(ix,4));
            assert(max(abs(pair.threshold(ix)-threshold))<1e-8);
        end
        % A changed future experimental input cannot change today's forecast.
        model=subject.models{source}; inputs=bse_inputs(sequence,model,source);
        if source==1, planned=bse_inputs(sequence,model,2); else, planned=inputs; end
        current=800; changed=inputs; changed(current+1:end,:)=.5*changed(current+1:end,:);
        a=nmpc_forecast(sequence,inputs,planned,current,pair.options);
        b=nmpc_forecast(sequence,changed,planned,current,pair.options);
        assert(isequaln(a,b));
    end
end
% Check rising-ramp and immediate-withdrawal semantics inside the Hill solver.
sequence=subject.sequences{1}; model=subject.models{1};
base=bse_model(model,subject.balanceController,sequence.standingAngle,true);
plant=instant_bse_actuator(base,model,subject.balanceController,sequence.standingAngle);
inputs=bse_inputs(sequence,model,1);
[state,~]=plant.initialize(inputs(1,:),sequence.mass(1),sequence.waiting(1),sequence.initialState);
forecast=struct('activation',inputs(1,:),'mass',sequence.mass(1),'waiting',sequence.waiting(1));
options=settings.controllers{1}; options.predictionDt=.01; options.predictionSubsteps=10;
[~,stepped]=plant.step(state,inputs(1,:),sequence.mass(1),sequence.waiting(1),2,.01,10);
predicted=plant.preview(state,50,forecast,options);
assert(stepped(25)==2 && max(abs(stepped-predicted),[],'all')<1e-7);
[~,originalRise]=base.step(state,inputs(1,:),sequence.mass(1),sequence.waiting(1),2,.01,10);
assert(max(abs(stepped-originalRise),[],'all')<1e-7,'Increasing torque no longer follows the original ramp.');
fallingState=state; fallingState(29)=50;
[~,withdrawn]=plant.step(fallingState,inputs(1,:),sequence.mass(1),sequence.waiting(1),0,.01,10);
predicted=plant.preview(fallingState,0,forecast,options);
[~,ramped]=base.step(fallingState,inputs(1,:),sequence.mass(1),sequence.waiting(1),0,.01,10);
assert(withdrawn(25)==0 && max(abs(withdrawn-predicted),[],'all')<1e-7);
assert(abs(withdrawn(2)-ramped(2))>1e-6,'Withdrawal still ramps inside the integrator.');
fprintf('PASS: clipped references, fixed tunings, shared 80-N-m caps, 200-N-m/s increase, instantaneous withdrawal, 0.5-N-m reserve and causal forecasts.\n');
disp(array2table(torqueExcess,'VariableNames',{'ExperimentalMaximumExcessNm','MAPMaximumExcessNm'}));
assert(all(torqueExcess<=1e-8,'all'),'Comparison:Torque','Some simulated muscle torques exceed the threshold with the requested 0.5-N-m reserve.');
disp(array2table(poseMaximum,'VariableNames',{'ExperimentalMaximumThetaRad','MAPMaximumThetaRad'}));
assert(all(poseMaximum<1.8,'all'),'Comparison:Pose', ...
    'Torque/reference checks pass, but recovery-angle sanity still fails. See the full unmodified angle plots.');
end
