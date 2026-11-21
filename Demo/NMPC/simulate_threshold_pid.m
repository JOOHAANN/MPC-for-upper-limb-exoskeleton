function run = simulate_threshold_pid(plant, sequence, inputs, threshold, options, pid)
%SIMULATE_THRESHOLD_PID Reactive, unilateral PID on measured back torque.
% The controller is armed only after the actual torque exceeds the threshold.
% It receives no MAP preview, future activation or unassisted trajectory.
% Velocity-form PID retains the reachable command as controller memory.
% Saturating that memory prevents windup without freezing the entire I step.
% The derivative acts on filtered measurement, avoiding threshold-change kicks.

%% Initialize the shared Hill model and the PID memory
count = size(inputs,1);
dt = sequence.dt;
stride = round(options.controlPeriod/dt);
period = stride*dt;
[state, first] = plant.initialize(inputs(1,:),sequence.mass(1), ...
    sequence.waiting(1),sequence.initialState);
run.actual = zeros(count,27);
run.actual(1,:) = first;
run.command = zeros(count,1);
run.armed = false(count,1);
run.integralNm = zeros(count,1);
run.solveSeconds = zeros(count,1);
integral = 0;
derivative = 0;
previousTorque = first(4);
previousError = 0;
previousDerivative = 0;
command = 0;
armed = false;

%% Read current torque, update PID, and advance the same continuous state
for k = 2:count
    if mod(k-2,stride) == 0
        tick = tic;
        measured = run.actual(k-1,4);
        previousCommand=command;
        error = measured-threshold(k-1);
        decay = exp(-period/pid.derivativeTime);
        derivative = decay*derivative+(1-decay)*(measured-previousTorque)/period;
        previousTorque = measured;
        if sequence.waiting(k-1)
            command = 0;
            integral = 0;
            derivative = 0;
            armed = false;
            previousError = 0;
            previousDerivative = 0;
        else
            if error > 0, armed = true; end
            if armed
                increment = pid.Kp*(error-previousError) ...
                    + pid.Ki*period*error ...
                    + pid.Kd*(derivative-previousDerivative);
                reachable = min(options.capNm,state(29)+options.rateNmPerSec*period);
                command = max(0,min(reachable,command+increment));
                % Equivalent I contribution after saturation tracking. It may
                % be negative; it is a controller state, not physical torque.
                integral = command-pid.Kp*error-pid.Kd*derivative;
                previousError = error;
                previousDerivative = derivative;
                if command < 1e-8 && error <= 0 && abs(integral) < 1e-6
                    armed = false;
                end
            else
                command = 0;
            end
        end
        % Optional command shaping smooths withdrawal, without delaying onset
        % or changing the actuator's instantaneous-decrease capability.
        if ~sequence.waiting(k-1) && isfield(pid,'releaseTime') && pid.releaseTime>0 ...
                && command<previousCommand
            command=command+(previousCommand-command)*exp(-period/pid.releaseTime);
            if command<.05,command=0;end
            integral=command-pid.Kp*error-pid.Kd*derivative;
        end
        run.solveSeconds(k) = toc(tick);
    end
    lower = 0; % Withdrawal is instantaneous; only increasing torque is ramped.
    upper = min(options.capNm,state(29)+options.rateNmPerSec*dt);
    applied = max(lower,min(upper,command));
    [state,run.actual(k,:)] = plant.step(state,inputs(k,:),sequence.mass(k), ...
        sequence.waiting(k),applied,dt,options.actualSubsteps);
    run.command(k) = command;
    run.armed(k) = armed;
    run.integralNm(k) = integral;
end
run.threshold = threshold;
run.pid = pid;
end
