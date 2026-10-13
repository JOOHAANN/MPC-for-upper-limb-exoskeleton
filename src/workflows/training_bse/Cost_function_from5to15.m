function Cost = Cost_function_from5to15(m1, theta_real_1, Torque_real_1, Smooth_activation_1, Smooth_activation_front, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, LMT0, front_a, front_b, front_c)

    dtheta_real_1 = [0; diff(theta_real_1) ./ Ts];
    ddtheta_real_1 = [0; diff(dtheta_real_1) ./ Ts];
    Tac = zeros(length(Smooth_activation_1), 1)';

    X1 = Torque_function_constant(m1, Tac, theta_real_1, theta_real_1, Smooth_activation_1, Smooth_activation_front, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, LMT0, front_a, front_b, front_c);

    % Y1 = Torque_function_constant(m1, 50, theta_real_1, Torque_real_1, Smooth_activation_1, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, LMT0);

    % X2 = Torque_function_constant(m2, 0, theta_real_2, Torque_real_2, Smooth_activation_2, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, Bias, LMT0);

    % Y2 = Torque_function_constant(m2, 80, theta_real_2, Smooth_activation_2, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, d, bia);

    % Zero_activation = zeros(length(Smooth_activation_1),1)';
    % Z1 = Torque_function_constant(m1, 0, theta_real_1, Zero_activation, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, Bias, LMT0);
    % Z2 = Torque_function_constant(m1, 80, theta_real_1, Zero_activation, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, LMT0);

    % Torque_real_50_2 = Torque_real_2 - 80;
    %
    % if Torque_real_50_2 < 5
    %     Torque_real_50_2 = 5;
    % end

    % Torque_mean = Z(:,1)/length(activation_zero);

    % x =  1.5*norm(X1(:,1) - Torque_real_1) + 1.5 * norm(Y1(:,1) - Torque_real_50_1) + norm(X2(:,1) - Torque_real_2) + 0.8 * norm(Y2(:,1) - Torque_real_50_2);
    % LM_LoM = mean(Z1(:,10));
    % if LM_LoM>1.35 || LM_LoM < 1.15 || any(isnan(LM_LoM)) || any(isinf(LM_LoM)) || ~isreal(LM_LoM)
    %     if LM_LoM > 1.35
    %     LM_LoM = LM_LoM + 1e3 * (LM_LoM - 1.35);  % 惩罚与偏差成比例
    %     elseif LM_LoM < 1.15
    %     LM_LoM = LM_LoM + 1e3 * (1.15 - LM_LoM);  % 惩罚与偏差成比例
    %     end
    % end
    % LM_LoM_P = LM_LoM;
    % if LM_LoM > 1.3
    %     LM_LoM_P = LM_LoM + 1e9 * (LM_LoM - 1.3);  % 惩罚与偏差成比例
    % elseif LM_LoM < 1.2
    %     LM_LoM_P = LM_LoM + 1e9 * (1.2 - LM_LoM);  % 惩罚与偏差成比例
    % elseif any(isnan(LM_LoM)) || any(isinf(LM_LoM)) || ~isreal(LM_LoM)
    %     LM_LoM_P = 10000;
    % end
    % Cost =  norm(X1(:,1) - Torque_real_1) + var(Z1(:,1)) + var(Z2(:,1));

    % Cost =  2.5*norm(X1(1:8000,1) - Torque_real_1(1:8000)) + 200*norm(X1(:,2) - theta_real_1) + 25000*var(X1(1:8000,2)) + 8*norm(X1(8109:13074,1) - Torque_real_1(8109:13074));
    % Cost_torque = 3.5*norm(X1(8109:13074,1) - Torque_real_1(8109:13074));
    % Cost_theta = 180*norm(X1(8109:13074,2) - theta_real_1(8109:13074));
    %
    % Cost_stable = 0.5*norm(X1(1:8000,1) - Torque_real_1(1:8000)) + 5000*var(X1(1:8000,2));
    % Cost_extra = 2.5*max(Y1(:,1));
    % Cost =  Cost_torque + Cost_theta + Cost_stable + Cost_extra;

    k = 1:length(Torque_real_1);
    valid_idx = (k > 7757 & k < 13853) | (k > 20956 & k < 27639);
    % valid_idx = (k > 10572 & k < 13853) | (k > 24037 & k < 27639) | (k > 39333 & k < 41734);
    % valid_idx = (k > 2767 & k < 8917) ;
    %| (k > 36908 & k < 41734);

    % valid_idx = (k > 7757 & k < 13853) ;

    Cost_torque = 2 * norm(X1(valid_idx, 1) - Torque_real_1(valid_idx));
    Cost_theta = 80 * norm(X1(valid_idx, 2) - theta_real_1(valid_idx)) + 80 * norm(X1(23500:24820, 2) - theta_real_1(23500:24820));
    Cost_dtheta = 0 * norm(X1(:, 3) - dtheta_real_1);
    Cost_ddtheta = 0 * norm(X1(:, 4) - ddtheta_real_1);

    % Cost_stable = 0*norm(X1(1:8000,1) - Torque_real_1(1:8000)) + 0*norm(X1(1:8000,2) - theta_real_1(1:8000)) + 0*var(X1(1:8000,2)) + 0*norm(Y1(13200:end,2) - theta_real_1(13200:end));
    Cost =  Cost_torque + Cost_theta + Cost_dtheta + Cost_ddtheta;

    % disp(['Current value of Cost_torque: ', num2str(Cost_torque)]);
    %% make the system stable
    % activation_zero = zeros(length(Smooth_activation),1)';
    % Y = Torque_function_constant(Real_theta, activation_zero, damping, alpha0, Ts, L1, R, LoM_nominal, LsT, FoM, L3, I);
    % x = 25*norm(Y(:,4) - 0) + norm(X(:,1) - Torque_real);
    %% make the system stable

    % if any(X(:,2) > 2.5) || any(X(:,2) < 0)
    %     x = 1e6;
    % else
    %     x = norm(X(:,1) - Torque_real) + 250*norm(X(:,2) - theta);
    % end

    if isnan(Cost)
        Cost = 1e9;
    end

    if isinf(Cost)
        Cost = 1e9;
    end

end
