% Inspect one MAP torque peak without rerunning an entire closed-loop trial.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','NMPC');addpath(demo,fullfile(root,'Demo','BSE'));
d=load(fullfile(demo,'output','no_human_PD_repair_priority','group1_source2.mat'));
subjects=bse_subjects('528'); subject=subjects{1}; q=d.sequence;
plant=bse_without_posture(subject.models{2},subject.balanceController,q.standingAngle,.8);
inputs=d.inputs; [state,~]=plant.initialize(inputs(1,:),q.mass(1),q.waiting(1),q.initialState);
solver=optimoptions('quadprog','Algorithm','interior-point-convex','Display','off', ...
    'ConstraintTolerance',1e-8,'OptimalityTolerance',1e-7,'MaxIterations',100);
rows=[];
for k=2:1001
    u=d.pair.nmpc.actual(k,25);
    [state,~]=plant.step(state,inputs(k,:),q.mass(k),q.waiting(k),u,q.dt,10);
    if ~ismember(k,[751 851 951 1001]), continue, end
    o=d.options; forecast=no_posture_forecast(q,inputs,inputs,k,o);
    query=min(numel(q.time),k+(1:o.points)'*o.predictionDt/q.dt);
    ref.theta=interp1((1:numel(q.time))',d.task.theta,query);
    ref.velocity=interp1((1:numel(q.time))',d.task.velocity,query);
    ref.threshold=interp1((1:numel(q.time))',d.pair.threshold,query);
    ref.torque=min(interp1((1:numel(q.time))',d.task.torque,query),ref.threshold);
    warm=state(29)*ones(o.blocks,1);
    for iteration=1:4
        [command,plan,flag,diagnostic]=no_posture_control_update(plant,state,warm,forecast,ref,o,solver);
        rows=[rows;k,iteration,command,flag,diagnostic]; %#ok<AGROW>
        warm=plan;
    end
end
disp(array2table(rows,'VariableNames',{'Index','Iteration','Command','Flag','PhaseFlag','MotionFlag', ...
    'TorqueSlackNm','AngleSlackRad','PredictedMaximumExcessNm','MaximumAngleErrorRad', ...
    'FirstTorqueNm','FirstThetaRad','Violation'}));
writematrix(rows,fullfile(folder,'optimizer_diagnostics.csv'));
