function test_pid_damping_528
%TEST_PID_DAMPING_528 Check the retuned reactive PID and unchanged NMPC traces.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'Demo','NMPC'); addpath(folder);
old=load(fullfile(root,'experiments','pid_damping_528','before_comparison_528.mat'));
new=load(fullfile(folder,'output','pid_nmpc_528','threshold_tracking','pid_nmpc_comparison_528.mat'));
settings=load(fullfile(folder,'data','threshold_tracking_settings_528.mat'));
assert(isequaln(settings.pid,new.results.pid));
assert(settings.pid.Ki<old.results.pid.Ki);
for group=1:3
    for source=1:2
        a=new.results.groups{group}.runs{source};
        b=old.results.groups{group}.runs{source};
        assert(isequaln(a.nmpc.actual,b.nmpc.actual),'NMPC traces changed.');
        assert(isequaln(a.baseline,b.baseline) && isequaln(a.threshold,b.threshold));
        assert(isequaln(a.pid.pid,settings.pid));
        assert(all(isfinite(a.pid.actual),'all'));
        u=a.pid.actual(:,25); dt=new.results.groups{group}.sequence.dt;
        assert(min(u)>=0 && max(u)<=80+1e-8);
        assert(max(diff(u))/dt<=200+1e-6);
        assert(max(abs(u(2:end)-min(a.pid.command(2:end),u(1:end-1)+200*dt)))<1e-8);
    end
end
% A positive error must not accumulate integral while the motor cannot follow.
plant=struct('initialize',@initialize_test,'step',@step_test);
sequence=struct('dt',.01,'mass',zeros(6,1),'waiting',[0;0;0;1;1;1],'initialState',[]);
options=struct('controlPeriod',.02,'capNm',80,'rateNmPerSec',200,'actualSubsteps',10);
pid=struct('Kp',1,'Ki',60,'Kd',0,'derivativeTime',.1);
run=simulate_threshold_pid(plant,sequence,zeros(6,2),zeros(6,1),options,pid);
assert(run.integralNm(2)==0,'The integral winds up against the rising ramp.');
assert(run.command(2)>0 && run.actual(2,25)==2);
assert(run.command(6)==0 && run.integralNm(6)==0,'Standing reset failed.');
fprintf('PASS: one PID across 18 cases, finite trajectories, shared actuator constraints, rise-rate anti-windup, standing reset and unchanged NMPC/baseline/threshold traces.\n');
end

function [state,out]=initialize_test(varargin)
state=zeros(1,29); out=zeros(1,27); out(4)=10;
end

function [state,out]=step_test(state,activation,mass,waiting,u,dt,steps) %#ok<INUSD>
state(29)=u; out=zeros(1,27); out(4)=10; out(25)=u;
end
