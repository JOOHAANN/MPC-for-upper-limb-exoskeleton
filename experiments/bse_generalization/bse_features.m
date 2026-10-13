function [lumped, physics, residual] = bse_features(trial)
%BSE_FEATURES Build fixed, target-free bases for torque estimation.
% The physics basis is a reduced estimator inspired by gravity, inertia,
% damping, and muscle activation. Its regression coefficients are not a
% physiological identification of the original full BSE dynamics.
    q = trial.q;
    v = trial.dq;
    acc = trial.ddq;
    e = trial.emg;
    back = mean(e(:, 3:8), 2);
    front = mean(e(:, 1:2), 2);
    m = trial.loadKg;
    lumped = back;
    physics = [back, front, sin(q), m * sin(q), ...
        v(:, 1), acc(:, 1), m * acc(:, 1)];
    % Keep channels separate to retain co-contraction and left/right effects.
    % Nominal load is context, not an assertion of continuous object contact.
    residual = [e, trial.emgPast, q, v, back .* q, back .* v, ...
        front .* q, e.^2, q.^2, v .* abs(v), ...
        m * e, m * q.^2, back.^2 .* q];
end
