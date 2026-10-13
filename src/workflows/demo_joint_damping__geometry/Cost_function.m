function Cost = Cost_function(theta_real_1, Torque_real_1, Smooth_activation_1, Ts, M, L1, L2, phi, Torque0, LMT0)
    % COST_FUNCTION Evaluate the calibration objective for this model variant.
    % Use prepare_workflow to select matching code and input data.

    Tac = zeros(length(Smooth_activation_1), 1)';

    X1 = Torque_function_constant(Tac, theta_real_1, Smooth_activation_1, Ts, M, L1, L2, phi, Torque0, LMT0);

    k = 1:length(Torque_real_1);
    valid_idx = (k > 7757 & k < 13853) | (k > 20956 & k < 27639) | (k > 34628 & k < 38485);

    Cost_torque = norm(X1(valid_idx, 1) - Torque_real_1(valid_idx));

    % Cost_stable = 0*norm(X1(1:8000,1) - Torque_real_1(1:8000)) + 0*norm(X1(1:8000,2) - theta_real_1(1:8000)) + 0*var(X1(1:8000,2)) + 0*norm(Y1(13200:end,2) - theta_real_1(13200:end));
    Cost =  Cost_torque;

    % disp(['Current value of Cost_torque: ', num2str(Cost_torque)]);

    %% make the system stable

    if isnan(Cost)
        Cost = 1e9;
    end

    if isinf(Cost)
        Cost = 1e9;
    end

end
