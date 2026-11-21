function test_pid_reachable_memory
%TEST_PID_REACHABLE_MEMORY Test anti-windup accounting, not human dynamics.
% A fixed torque sensor isolates the I/rate-limit interaction. It is not a
% musculoskeletal step-response or settling-time performance experiment.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'Demo','NMPC'));
q=struct('dt',.01,'mass',10*ones(201,1),'waiting',zeros(201,1),'initialState',[1.5 0]);
plant=struct('initialize',@initialize,'step',@advance);
o=struct('controlPeriod',.02,'capNm',80,'rateNmPerSec',200,'actualSubsteps',10);
pid=struct('Kp',.413254656908,'Ki',6.35202882737,'Kd',0,'derivativeTime',.1,'releaseTime',0);
r=simulate_threshold_pid(plant,q,50*ones(201,2),100*ones(201,1),o,pid);
u=r.actual(:,25);
assert(max(u)<=80+1e-10 && max(diff(u))/.01<=200+1e-8);
assert(find(u>=80-1e-9,1)*.01<=.5 && u(end)>=80-1e-9, ...
    'The integral remains stuck below the amplitude cap because of the rise limit.');
% Without saturation, velocity and position PID forms must be equivalent.
pid.Kp=.5;pid.Ki=.1;
r=simulate_threshold_pid(plant,q,5*ones(201,2),100*ones(201,1),o,pid);
assert(abs(r.actual(end,25)-3.5)<1e-9);
fprintf('PASS: reachable-output PID builds support through a rate limit and matches unsaturated PI algebra.\n');
end

function [state,out]=initialize(input,~,~,initial)
state=zeros(69,1);state(1:2)=initial;out=zeros(1,27);
out(1)=initial(1);out(4)=100+input(1);out(3)=out(4);
end

function [state,out]=advance(state,input,~,~,applied,~,~)
state(29)=applied;out=zeros(1,27);out(1)=state(1);
out(4)=100+input(1);out(3)=out(4);out(25)=applied;
end
