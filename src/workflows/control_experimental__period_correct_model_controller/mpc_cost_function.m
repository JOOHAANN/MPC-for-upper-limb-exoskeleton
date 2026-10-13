function J = mpc_cost_function(z, InitialConditions, N, Q, R, r, Ts, damping, alpha0, activation, coeffs)
    J = 0;
    u_Tac = z;
    x = nonlinear_system_period(u_Tac, InitialConditions, activation, damping, alpha0, Ts, coeffs);
    for i = 1:N
        J = J + (x(i, 1) - r)' * Q * (x(i, 1) - r) + R * u_Tac(i)^2;
    end
    % J = real(J);
end
