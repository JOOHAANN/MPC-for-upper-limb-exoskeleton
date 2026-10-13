function J = mpc_cost_function_increment(m, u0, du, n, InitialConditions, PH, CH, Q, R, W, r, damping, alpha0, Ts, activation, theta_real, dtheta_real, coeffs)
    % MPC_COST_FUNCTION_INCREMENT Evaluate the incremental nonlinear predictive-control objective.
    % Use prepare_workflow to select matching code and input data.
    J = 0;
    u = u0;
    u_Tac = zeros(PH, 1);
    CHn = CH * n;
    du_full = repelem(du, n);
    for i = 1:PH
        if i <= CHn
            u = u + du_full(i);
            u_Tac(i) = u;
        else
            u_Tac(i) = u;
        end
    end
    x = nonlinear_system_period(m, u_Tac, InitialConditions, activation, theta_real, dtheta_real, damping, alpha0, Ts, coeffs);
    for i = 1:PH
        error = x(i, 1) - r(i);
        if error > 0
            penalty = (10 * error)^2;
        else
            penalty = 0;
        end
        if i > CHn
            J = J + penalty' * Q * penalty + R * u_Tac(i)^2;
        elseif i < 3
            J = 0;
        else
            J = J + penalty' * Q * penalty + R * u_Tac(i)^2 + W * du_full(i)^2;
        end
    end
end
