function Cost = Cost_function(theta_real_1, Torque_real_1, Smooth_activation_1, Smooth_activation_front, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, zeta, LMT0)
    % COST_FUNCTION Evaluate the calibration objective for this model variant.
    % Use prepare_workflow to select matching code and input data.

    dtheta_real_1 = [0; diff(theta_real_1) ./ Ts];
    ddtheta_real_1 = [0; diff(dtheta_real_1) ./ Ts];
    Tac = zeros(length(Smooth_activation_1), 1)';

    X1 = Torque_function_constant(Tac, theta_real_1, Torque_real_1, Smooth_activation_1, Smooth_activation_front, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, zeta, LMT0);

    k = 1:length(Torque_real_1);
    % valid_idx = (k > 7757 & k < 13853) | (k > 20956 & k < 27639) | (k > 36895 & k < 41683);
    valid_idx_354 = (k > 5321 & k < 11740) | (k > 22037 & k < 26059) | (k > 36732 & k < 41109);
    valid_idx_435 = (k > 52761 & k < 58381) | (k > 66452 & k < 72417) | (k > 84099 & k < 87884);
    valid_idx_543 = (k > 98004 & k < 102357) | (k > 113394 & k < 118437) | (k > 127817 & k < 133717);

    valid_idx = valid_idx_354 | valid_idx_435 | valid_idx_543;

    Cost_torque = 1 * norm(X1(valid_idx, 1) - Torque_real_1(valid_idx));

    % Cost_theta = calculate_cost_threshold(theta_real_1(valid_idx), X1(valid_idx,2));

    Cost_theta = calculate_cost_hybrid(theta_real_1(valid_idx), X1(valid_idx, 2));

    % Cost_stable = 0*norm(X1(1:8000,1) - Torque_real_1(1:8000)) + 0*norm(X1(1:8000,2) - theta_real_1(1:8000)) + 0*var(X1(1:8000,2)) + 0*norm(Y1(13200:end,2) - theta_real_1(13200:end));
    Cost =  Cost_torque + Cost_theta;

    % disp(['Current value of Cost_torque: ', num2str(Cost_torque)]);

    %% make the system stable

    if isnan(Cost)
        Cost = 1e9;
    end

    if isinf(Cost)
        Cost = 1e9;
    end

end
