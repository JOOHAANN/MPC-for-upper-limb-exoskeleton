function run = simulate_threshold_nmpc(plant, sequence, inputs, planned, threshold, options)
%SIMULATE_THRESHOLD_NMPC Predictive torque tracking with smooth assistance release.
% Two small QPs first minimize attainable torque slack, then track the entire
% upper-clipped unassisted torque while preserving the predicted motion.
% No measured reference trajectory enters the controller.
% A source-specific release setting prevents rapid loss of support between
% action peaks. Both controllers retain the same physical actuator limits.
assert(strcmp(options.inputForecast,'arx') && ~options.mapPhaseObserver, ...
    'NMPC:Forecast', 'Use the trained MAP-plus-ARX forecast configuration.');
solver.motion = optimoptions('quadprog','Algorithm','active-set','Display','off', ...
    'OptimalityTolerance',1e-9,'ConstraintTolerance',1e-9,'MaxIterations',120);
% The slack QP is nearly semidefinite; interior point avoids active-set cycling.
solver.feasibility = optimoptions('quadprog','Algorithm','interior-point-convex', ...
    'Display','off','OptimalityTolerance',1e-9,'ConstraintTolerance',1e-9,'MaxIterations',80);
%% Preallocate the simulation arrays and initialize the complete Hill state
count = size(inputs,1);
dt = sequence.dt;
[state, first] = plant.initialize(inputs(1,:),sequence.mass(1), ...
    sequence.waiting(1),sequence.initialState);
nominalState = state;
run.actual = zeros(count,27);
run.nominal = zeros(count,27);
run.actual(1,:) = first;
run.nominal(1,:) = first;
run.target = zeros(count,1);
run.target(1) = min(first(4),threshold(1));
run.solverExitFlag = nan(count,1);
run.solveSeconds = zeros(count,1);
run.command = zeros(count,1);
warm = zeros(options.blocks,1);
command = 0;
stride = round(options.controlPeriod/dt);
nextProgress = .2;
%% Receding-horizon simulation: forecast, optimize, apply and advance
for k = 2:count
    if mod(k-2,stride) == 0
        tick = tic;
        current = k-1;
        forecast = nmpc_forecast(sequence,inputs,planned,current,options);
        online=options;
        if sequence.waiting(current), online.releaseRateNmPerSec=options.rateNmPerSec; end
        % Actions and standing waits use the same clipped-reference objective.
        identicalWait=isfield(options,'currentCalibratedModel') && options.currentCalibratedModel ...
            && all(forecast.waiting>.5) && state(29)==0 ...
            && max(abs(state(1:28)-nominalState(1:28)))<1e-10;
        if identicalWait
            command=0;flag=1;plan=zeros(options.blocks,1);
        else
            [command, flag, plan] = control_update(plant,state,nominalState,warm, ...
                threshold(current),forecast,online,solver);
        end
        % The governor changes the command, never the simulated joint angle.
        % Inspect the optimized future release, rather than holding today's
        % command constant throughout the entire prediction horizon.
        plan(1)=command;
        if options.angleGovernorUpperRad > 0 && (command > 1e-6 || state(29) > 1e-6)
            % Recovery posture protection acts near upright. It must not
            % withdraw needed support during the deep-flexion torque peaks.
            hasMotionLimit=isfield(options,'maximumAngleLeadRad') && options.maximumAngleLeadRad>0 ...
                && state(1)>sequence.standingAngle-.15 && state(2)>0;
            if hasMotionLimit
                nominalPreview=controller_preview(plant,nominalState,zeros(options.blocks,1),forecast,options);
            end
            for attempt = 1:8
                predicted = controller_preview(plant,state,plan,forecast,options);
                followsMotion=true;
                if hasMotionLimit
                    followsMotion=~isempty(nominalPreview) && ~isempty(predicted) ...
                        && max(predicted(:,1)-nominalPreview(:,1))<=options.maximumAngleLeadRad;
                end
                if ~isempty(predicted) && max(predicted(:,1)) <= options.angleGovernorUpperRad && followsMotion
                    break
                end
                if command <= 1e-6, break, end
                if command < .1, plan(:)=0; else, plan=.5*plan; end
                command=plan(1);
            end
        end
        % Recede the optimized sequence instead of replacing it by a flat hold.
        % A fractional shift retains the planned onset and withdrawal timing.
        controlPoints=options.points;
        if isfield(options,'controlPoints'), controlPoints=options.controlPoints; end
        blockTime=options.predictionDt*controlPoints/options.blocks;
        grid=(0:options.blocks-1)'*blockTime;
        warm=interp1(grid,plan,grid+options.controlPeriod,'linear',plan(end));
        warm=max(0,min(options.capNm,warm));
        run.solverExitFlag(k) = flag;
        run.solveSeconds(k) = toc(tick);
    end
    lower = 0; % No torque-decrease rate limit.
    upper = min(options.capNm,state(29)+options.rateNmPerSec*dt);
    applied = max(lower,min(upper,command));
    run.command(k)=command;
    [state, run.actual(k,:)] = plant.step(state,inputs(k,:),sequence.mass(k), ...
        sequence.waiting(k),applied,dt,options.actualSubsteps);
    [nominalState, run.nominal(k,:)] = plant.step(nominalState,inputs(k,:),sequence.mass(k), ...
        sequence.waiting(k),0,dt,options.actualSubsteps);
    run.target(k) = min(run.nominal(k,4),threshold(k));
    if k/count >= nextProgress
        fprintf('  %3.0f%% of continuous sequence\n',100*k/count);
        nextProgress = nextProgress+.2;
    end
end
run.threshold = threshold;
run.inputs = inputs;
end

function [command, flag, plan] = control_update(plant, state, nominalState, warm, threshold, forecast, options, solver)
%% Set the actuator bounds and prediction uncertainty reserve
n = numel(warm);
plan=warm;
current = state(29);
cap = options.capNm;
priority = options.torqueConstraintPriority;
referenceLimit = max(0,threshold);
safetyLimit = max(0,referenceLimit-options.backTorqueMarginNm);
nominal = controller_preview(plant,nominalState,zeros(n,1),forecast,options);
limit = options.rateNmPerSec*options.controlPeriod;
releaseRate=options.releaseRateNmPerSec;
% Withdraw faster when the current unassisted shadow model needs support,
% yet assisted muscle torque is already well below its limit. Only current
% model states enter this rule; no measured reference or future EMG is used.
if isfield(options,'releaseTrackingBandNm') && options.releaseTrackingBandNm>0
    backCapacity=plant.physics.back.capacity;
    assistedBack=backCapacity*state(16);
    unassistedBack=backCapacity*nominalState(16);
    currentTarget=min(unassistedBack,referenceLimit);
    if assistedBack<currentTarget-options.releaseTrackingBandNm
        releaseRate=options.rateNmPerSec;
    end
end
releaseLimit=releaseRate*options.controlPeriod;
if options.instantRelease
    releaseLimit=cap;
end
if isempty(nominal)
    command = max(0,current-limit); flag = -3; return
end
guard = options.predictionErrorGuardNm(:);
if ~isempty(options.predictionErrorGuardRelative)
    guard = 2+options.predictionErrorGuardRelative(:)*referenceLimit/.75;
end
floorNm = options.predictionErrorGuardFloorNm;
guard = floorNm+options.predictionErrorGuardScale*max(guard-floorNm,0);
if options.predictionErrorGuardCapNm > 0
    guard = min(guard,options.predictionErrorGuardCapNm);
end
ceiling = safetyLimit-guard;
constraintPoints = min(options.points,max(1,round( ...
    options.torqueConstraintHorizonSeconds/options.predictionDt)));
if max(nominal(:,4)-ceiling) < -5 && current < 1e-8 ...
        && (~priority || max(abs(state(1:2)-nominalState(1:2))) < 1e-7)
    command = 0; flag = 1; return
end
lower = zeros(n,1); upper = cap*ones(n,1);
lower(1) = max(0,current-releaseLimit); upper(1) = min(cap,current+limit);
controlPoints=options.points;
if isfield(options,'controlPoints'), controlPoints=options.controlPoints; end
blockLimit = options.rateNmPerSec*options.predictionDt*controlPoints/n;
blockReleaseLimit=releaseRate*options.predictionDt*controlPoints/n;
if options.instantRelease
    blockReleaseLimit=cap;
end
warm = max(0,min(cap,warm));
warm(1) = max(lower(1),min(upper(1),warm(1)));
for j = 2:n
    warm(j) = max(max(0,warm(j-1)-blockReleaseLimit),min(min(cap,warm(j-1)+blockLimit),warm(j)));
end
%% Linearize the nonlinear Hill prediction around the previous command
base = controller_preview(plant,state,warm,forecast,options);
if isempty(base), command = max(0,current-limit); flag = -3; return, end
r0 = objective_residual(base,warm,current,nominal,referenceLimit,options);
jacobian = zeros(numel(r0),n);
torqueJacobian = zeros(options.points,n);
for j = 1:n
    trial = warm;
    delta = .5;
    if trial(j) >= cap-.5, delta = -.5; end
    trial(j) = trial(j)+delta;
    predicted = controller_preview(plant,state,trial,forecast,options);
    if ~isempty(predicted)
        torqueJacobian(:,j) = (predicted(:,4)-base(:,4))/delta;
        jacobian(:,j) = (objective_residual(predicted,trial,current,nominal,referenceLimit,options)-r0)/delta;
    end
end
hessian = jacobian'*jacobian+1e-8*eye(n);
gradient = jacobian'*r0;
difference = diff(eye(n),1,1);
matrix = [difference; -difference];
rhs = [blockLimit*ones(n-1,1);blockReleaseLimit*ones(n-1,1)];
initial = warm;
%% Solve the torque-feasibility problem before the motion objective
if priority
    torqueJacobian = torqueJacobian(1:constraintPoints,:);
    torqueRhs = (ceiling(1:constraintPoints)-base(1:constraintPoints,4)+torqueJacobian*warm)/30;
    phaseMatrix = [matrix,zeros(size(matrix,1),1); torqueJacobian/30,-ones(constraintPoints,1)];
    phaseRhs = [rhs;torqueRhs];
    [phase,~,phaseFlag] = quadprog(diag([2e-10*ones(n,1);2]),zeros(n+1,1), ...
        phaseMatrix,phaseRhs,[],[],[lower;0],[upper;inf],[],solver.feasibility);
    if phaseFlag <= 0 || any(~isfinite(phase))
        command = warm(1); flag = -2; return
    end
    slack = max(0,phase(end))+1e-5;
    matrix = [matrix;torqueJacobian/30];
    rhs = [rhs;torqueRhs+slack];
    initial = phase(1:n);
end
[solution,~,exitFlag] = quadprog(hessian,gradient-hessian*warm,matrix,rhs, ...
    [],[],lower,upper,initial,solver.motion);
feasible = exitFlag > 0 && ~isempty(solution) && all(isfinite(solution)) ...
    && all(solution >= lower-1e-7) && all(solution <= upper+1e-7) ...
    && all(difference*solution <= blockLimit+1e-7) ...
    && all(-difference*solution <= blockReleaseLimit+1e-7);
if priority && feasible
    feasible = all(torqueJacobian*solution <= 30*(torqueRhs+slack)+1e-5);
end
if ~feasible
    if priority
        solution = initial;
    else
        command = warm(1); flag = -2; return
    end
end
%% Check the proposed move against the nonlinear prediction
best = warm;
value = r0'*r0;
violation = max(max(base(1:constraintPoints,4)-ceiling(1:constraintPoints),0));
fractions = [1 .5];
if priority, fractions = [1 .5 .25 .125 .0625]; end
for fraction = fractions
    candidate = warm+fraction*(solution-warm);
    predicted = controller_preview(plant,state,candidate,forecast,options);
    if isempty(predicted), continue, end
    residual = objective_residual(predicted,candidate,current,nominal,referenceLimit,options);
    excess = max(max(predicted(1:constraintPoints,4)-ceiling(1:constraintPoints),0));
    if priority
        permitted=max(violation,options.nonlinearTorqueToleranceNm);
        better = excess < violation-1e-3 || (excess <= permitted+1e-3 && residual'*residual < value);
    else
        better = residual'*residual < value;
    end
    if better
        best = candidate; value = residual'*residual; violation = excess;
    end
end
command = max(lower(1),min(upper(1),best(1)));
plan=best;
flag = 1;
if ~feasible, flag = -2; end
end

function output=controller_preview(plant,state,controls,forecast,options)
% A zero-assistance tail exposes delayed recovery after the optimized moves.
% The same Hill equations simulate both the controlled and coasting portions.
if isfield(options,'controlPoints') && options.controlPoints<options.points
    pointsPerBlock=options.controlPoints/numel(controls);
    tailBlocks=(options.points-options.controlPoints)/pointsPerBlock;
    controls=[controls;zeros(tailBlocks,1)];
end
output=plant.preview(state,controls,forecast,options);
end

function residual = objective_residual(predicted, controls, current, nominal, threshold, options)
% Follow the complete upper-clipped unassisted torque, including its valleys.
% Safety reserves constrain feasibility; they do not change this reference.
target = min(nominal(:,4),threshold);
abdomenPenalty=max(predicted(:,5)-12,0)/15;
if isfield(options,'abdomenReferenceWeight')
    abdomenPenalty=sqrt(options.abdomenReferenceWeight)* ...
        max(predicted(:,5)-nominal(:,5)-options.abdomenAllowanceNm,0)/15;
end
residual = [sqrt(options.angleWeight)*(predicted(:,1)-nominal(:,1))/.10; ...
    sqrt(options.backTorqueWeight)*(predicted(:,4)-target)/30; ...
    abdomenPenalty; ...
    sqrt(options.assistanceEffortWeight)*controls/80; ...
    .20*diff([current;controls])/20];
if options.velocityWeight > 0
    residual = [residual; sqrt(options.velocityWeight)*(predicted(:,2)-nominal(:,2))];
end
end
