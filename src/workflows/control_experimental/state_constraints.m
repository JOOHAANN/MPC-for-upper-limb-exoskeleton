function [c, ceq] = state_constraints(z, InitialConditions, N, Ts, LMT_min, activation, coeffs, k_t)
    L2 = coeffs(1);
    L3 = coeffs(2);
    R = coeffs(3);
    c = [];
    LMT = zeros(N, 1);  % 预定义LMT数组的大小
    c = zeros(N, 1);  % 预定义c数组的大小
    for i = 1:N
        u_Tac = z(i);
        [MTF_next, theta_next, dtheta_next] = nonlinear_system(k_t, InitialConditions, activation(i), u_Tac, Ts, coeffs);
        LMT(i) = L3^2 + L2^2 + 2 * L2 * L3 * cos(theta_next) - R^2;
        c = [c
             LMT_min - LMT(i)
            ];
    end
    ceq = [];
end
