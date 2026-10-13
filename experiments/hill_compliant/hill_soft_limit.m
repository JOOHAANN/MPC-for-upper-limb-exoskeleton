function torque = hill_soft_limit(q, par)
    % HILL_SOFT_LIMIT Conservative elastic recovery outside a learned dead zone.
    % No target/reference angle is used. It does not clip or replace the state.
    % The 0.05-rad onset makes the spring torque and its slope continuous.
    lower = max(par(31) - q, 0);
    upper = max(q - par(32), 0);
    lower = lower.^2 ./ (lower + .05);
    upper = upper.^2 ./ (upper + .05);
    torque = par(30) * (tanh(par(29) * lower / par(30)) - tanh(par(29) * upper / par(30)));
end
