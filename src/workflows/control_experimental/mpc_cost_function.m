function J = mpc_cost_function(u0, du, InitialConditions, N, Q, R, r, damping, alpha0, Ts, activation, Theta_Tac0, dTheta_Tac0, coeffs)
    J = 0;
    u_Tac = u0 + du;
    x = nonlinear_system_period(u_Tac, InitialConditions, activation, Theta_Tac0, dTheta_Tac0, damping, alpha0, Ts, coeffs);
    for i = 1:N
        penalty = (x(i, 1) - r(i))^2;
        J = J + penalty' * Q * penalty + R * u_Tac(i)^2;
    end
end
