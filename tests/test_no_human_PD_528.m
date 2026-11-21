function test_no_human_PD_528
%TEST_NO_HUMAN_PD_528 Model-ablation, reference and matched-hardware contracts.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'Demo','NMPC'); addpath(folder,fullfile(root,'Demo','BSE'));
data=load(fullfile(folder,'output','no_human_PD_528','comparison_no_human_PD_528.mat'));
r=data.results; subjects=bse_subjects('528'); subject=subjects{1};
reference=load(fullfile(folder,'data','no_human_PD_task_528.mat'));
original=load(fullfile(folder,'output','pid_nmpc_528','threshold_tracking','pid_nmpc_comparison_528.mat'));
motor=load(fullfile(folder,'data','no_human_PD_motor_programs_528.mat'));
assert(height(r.metrics)==54,'Expected all 9 cases x 2 sources x 3 controllers.');
assert(r.protocol.humanPostureFeedback && ~r.protocol.humanMovementFeedback ...
    && r.protocol.standingBalanceFeedback && r.protocol.bidirectionalActuator);
assert(isequaln(r.protocol.standingParameters,subject.balanceController));
assert(~r.protocol.effortMatched && r.protocol.lowerTorqueNm==-80 && r.protocol.upperTorqueNm==80);
for group=1:3
    q=r.groups{group}.sequence;
    for source=1:2
        p=r.groups{group}.runs{source}; o=r.controllers{source};
        assert(isequaln(p.options,o)); assert(isequaln(p.pid.pid,r.pid));
        assert(isequaln(p.task,r.groups{group}.runs{3-source}.task));
        assert(isequaln(p.task,reference.tasks{group}));
        assert(isequaln(p.nmpc.inputs,motor.programs{group,source}));
        assert(isequaln(p.nmpc.planned,motor.plannedPrograms{group,source}));
        model=subject.models{source};
        assert(isequaln(motor.calibrationParameters{source},[model.parameters(:);model.neuralParameters(:)]));
        assert(max(abs(p.threshold-original.results.groups{group}.runs{source}.threshold))<1e-6);
        assert(max(abs(p.nmpc.target-min(p.task.torque,p.threshold)))<1e-10);
        if source==2
            ix=find(isfinite(p.nmpc.optimizerDiagnostics(:,7)));
            assert(max(abs(p.nmpc.optimizerDiagnostics(ix,7)-p.nmpc.actual(ix,4)))<1e-6);
            assert(max(abs(p.nmpc.optimizerDiagnostics(ix,8)-p.nmpc.actual(ix,1)))<1e-8);
        end
        for a={p.baseline,p.pid.actual,p.nmpc.actual}
            a=a{1}; assert(all(isfinite(a),'all'));
            assert(all(abs(a(:,22))<=a(:,20)*subject.balanceController(4)+1e-8), ...
                'Tracking correction is not limited to standing allocation.');
            assert(max(abs(a(:,3)-a(:,4)+a(:,5)))<1e-6);
        end
        for controlled={p.pid,p.nmpc}
            c=controlled{1}; u=c.actual(:,25); assert(max(abs(u))<=80+1e-8);
            moving=c.actual(:,20)<1e-10;
            if any(moving)
                assert(max(c.actual(moving,18:19)-c.inputs(moving,:),[],'all')<1e-8);
            end
            for k=2:numel(u)
                previous=u(k-1); requested=c.command(k);
                if previous*requested<0, previous=0; end
                expected=sign(requested)*min(abs(requested),abs(previous)+200*q.dt);
                assert(abs(u(k)-expected)<1e-8,'Actuator magnitude-rise limit violated.');
            end
        end
        for stage=1:3
            seg=q.segments{stage}; ix=seg.startIndex:seg.endIndex;
            assert(max(abs(p.threshold(ix)-reference.tasks{group}.thresholds(ix,source)))<1e-10);
        end
    end
end
% Predicting a torque reversal must match the same sequence of actual steps.
q=subject.sequences{1}; model=subject.models{1};
plant=bse_without_posture(model,subject.balanceController,q.standingAngle);
assert(all(plant.neuralParameters(4:5)==0) && plant.neuralParameters(7)==0 ...
    && isequal(plant.standingParameters(2:end),subject.balanceController(:)'));
inputs=motor.programs{1,1};
[state,~]=plant.initialize(inputs(1,:),q.mass(1),q.waiting(1),q.initialState);
state(29)=20;
forecast=struct('activation',repmat(inputs(1,:),3,1),'mass',q.mass(1)*ones(3,1), ...
    'waiting',q.waiting(1)*ones(3,1));
options=r.controllers{1}; options.predictionDt=.025; options.predictionSubsteps=25;
controls=[22;-5;-3]; predicted=plant.preview(state,controls,forecast,options);
actual=zeros(3,27);
for k=1:3
    previous=state(29); if previous*controls(k)<0, previous=0; end
    u=sign(controls(k))*min(abs(controls(k)),abs(previous)+5);
    [state,actual(k,:)]=plant.step(state,inputs(1,:),q.mass(1),q.waiting(1),u,.025,25);
end
assert(max(abs(actual-predicted),[],'all')<1e-7);
assert(isequal(actual(:,25),[22;-5;-3]));
% During pure movement, position/velocity error cannot change commands.
a=inputs(800,:);
[s,~]=plant.initialize(a,5,0,[1.0 0]);
changed=s; changed(1)=changed(1)+.03; changed(2)=changed(2)+.1;
[~,first]=plant.step(s,a,5,0,0,.01,10);
[~,second]=plant.step(changed,a,5,0,0,.01,10);
assert(max(abs(first(18:19)-a))<1e-12 && max(abs(second(18:19)-a))<1e-12);
[~,braked]=plant.step(s,a,5,0,-2,.01,10);
[~,supported]=plant.step(s,a,5,0,2,.01,10);
assert(abs(braked(18)-a(1))<1e-12 && braked(19)<=a(2)+1e-12);
assert(abs(supported(19)-a(2))<1e-12 && supported(18)<=a(1)+1e-12);
% Waiting must reproduce the original standing PD, including its saturation.
originalPlant=bse_model(model,subject.balanceController,q.standingAngle,true);
standingInput=repmat([.1 .02],101,1);waiting=ones(101,1);mass=zeros(101,1);
originalStanding=originalPlant.simulate(standingInput,mass,waiting,[1.35 .05],.01,10);
restoredStanding=plant.simulate(standingInput,mass,waiting,[1.35 .05],.01,10);
physical=[1:23 25:27];
assert(max(abs(originalStanding(:,physical)-restoredStanding(:,physical)),[],'all')<1e-7);
assert(max(abs(restoredStanding(:,22)))>1);
% Standing feedforward drive is preserved, and future observations are absent.
options=r.controllers{1}; planned=motor.plannedPrograms{1,1}; observed=motor.programs{1,1};
f=no_posture_forecast(q,observed,planned,1,options);
assert(any(f.activation(f.waiting>.99,:)>0,'all'));
altered=observed; altered(801:end,:)=.5*altered(801:end,:);
f1=no_posture_forecast(q,observed,planned,800,options);
f2=no_posture_forecast(q,altered,planned,800,options); assert(isequaln(f1,f2));
fprintf('PASS: movement PD disabled; original standing balance and transition retained, identical external PID/task across cases, shared signed 80-N-m caps and 200-N-m/s magnitude rise, immediate withdrawal, torque accounting and frozen thresholds.\n');
fprintf('PASS: no movement posture-to-command feedback, original standing disturbance response, causal command-space forecast and exact MAP one-step prediction/reversal parity.\n');
end
