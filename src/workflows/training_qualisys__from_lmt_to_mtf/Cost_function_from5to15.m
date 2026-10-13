function [Cost, LM_LoM] = Cost_function_from5to15(m1, m2, Torque_real_1, Torque_real_2, theta_real_1, theta_real_2, Smooth_activation_1, Smooth_activation_2, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, LMT0)

    X1 = Torque_function_constant(m1, 0, theta_real_1, Smooth_activation_1, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, LMT0);

    % Y1 = Torque_function_constant(m1, 80, theta_real_1, Smooth_activation_1, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, d, bia, LMT0);

    % X2 = Torque_function_constant(m2, 0, theta_real_2, Smooth_activation_2, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, d, bia, LMT0);
    %
    % Y2 = Torque_function_constant(m2, 80, theta_real_2, Smooth_activation_2, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, d, bia);

    % Zero_activation = zeros(length(Smooth_activation_2),1)';
    % Z = Torque_function_constant(m1, 0, theta_real_2, Zero_activation, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, LMT0);
    %
    % Torque_real_50_1 = Torque_real_1 - 80;
    %
    % if Torque_real_50_1 < 5
    %     Torque_real_50_1 = 5;
    % end

    % Torque_real_50_2 = Torque_real_2 - 80;
    %
    % if Torque_real_50_2 < 5
    %     Torque_real_50_2 = 5;
    % end

    % Torque_mean = Z(:,1)/length(activation_zero);

    % x =  1.5*norm(X1(:,1) - Torque_real_1) + 1.5 * norm(Y1(:,1) - Torque_real_50_1) + norm(X2(:,1) - Torque_real_2) + 0.8 * norm(Y2(:,1) - Torque_real_50_2);

    % if LM_LoM>1.35 || LM_LoM < 1.15 || any(isnan(LM_LoM)) || any(isinf(LM_LoM)) || ~isreal(LM_LoM)
    %     if LM_LoM > 1.35
    %     LM_LoM = LM_LoM + 1e3 * (LM_LoM - 1.35);  % 惩罚与偏差成比例
    %     elseif LM_LoM < 1.15
    %     LM_LoM = LM_LoM + 1e3 * (1.15 - LM_LoM);  % 惩罚与偏差成比例
    %     end
    % end
    % LM_LoM_P = LM_LoM;
    % if LM_LoM > 1.3
    %     LM_LoM_P = LM_LoM + 1e6 * (LM_LoM - 1.3);  % 惩罚与偏差成比例
    % elseif LM_LoM < 1.2
    %     LM_LoM_P = LM_LoM + 1e6 * (1.2 - LM_LoM);  % 惩罚与偏差成比例
    % elseif any(isnan(LM_LoM)) || any(isinf(LM_LoM)) || ~isreal(LM_LoM)
    %     LM_LoM_P = 10000;
    % end
    Cost =  norm(X1(:, 1) - Torque_real_1);

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
