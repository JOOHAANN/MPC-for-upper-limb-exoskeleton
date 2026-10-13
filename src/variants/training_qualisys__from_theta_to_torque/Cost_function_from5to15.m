function Cost = Cost_function_from5to15(m1, Torque_real_1, theta_real_1, Ts, M, L1, L2, phi)

    dtheta_real_1 = [0; diff(theta_real_1) ./ Ts];
    ddtheta_real_1 = [0; diff(dtheta_real_1) ./ Ts];

    X1 = Torque_function_constant(m1, Torque_real_1, theta_real_1, Ts, M, L1, L2, phi);

    % Y1 = Torque_function_constant(m1, 50, theta_real_1, Torque_real_1, Smooth_activation_1, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, LMT0);

    % Torque_mean = Z(:,1)/length(activation_zero);

    Cost_torque = norm(X1(:, 1) - Torque_real_1);
    % Cost_theta = 100*norm(X1(:,2) - theta_real_1);
    % Cost_dtheta = 10*norm(X1(:,3) - dtheta_real_1);
    % Cost_ddtheta = 100*norm(X1(:,4) - ddtheta_real_1);

    Cost =  Cost_torque;

    if isnan(Cost)
        Cost = 1e9;
    end

    if isinf(Cost)
        Cost = 1e9;
    end

end
