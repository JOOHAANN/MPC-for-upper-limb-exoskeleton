function [c, ceq] = nonlincon(m1, m2, Torque_real_1, Torque_real_2, theta_real_1, theta_real_2, Smooth_activation_1, Smooth_activation_2, damping, alpha0, Ts, coeffs, LMT0)
    % 调用 Cost_function_from5to15，并获取成本和 LM_LoM
    [~, LM_LoM] = Cost_function_from5to15(m1, m2, Torque_real_1, Torque_real_2, theta_real_1, theta_real_2, Smooth_activation_1, Smooth_activation_2, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), LMT0);

    % 施加不等式约束 (LM_LoM 需在 1.2 到 1.35 之间)

    c = [mean(LM_LoM) - 1.2; 1.35 - mean(LM_LoM)];

    % % 确保 c 为实数
    % if any(isnan(c)) || any(isinf(c)) || ~isreal(c)
    %     c = [-0.2; 0.2];
    % else
    %     c = [mean(LM_LoM) - 1.2; 1.35 - mean(LM_LoM)];
    % end

    % 没有等式约束
    ceq = [];

    % 确保 ceq 为实数
    if any(isnan(ceq)) || any(isinf(ceq)) || ~isreal(ceq)
        error('约束函数的等式返回了无效的值 (NaN, Inf, 或复数)');
    end
end
