function [LMT, MTF, Torque, MA] = Control_model(activation, theta, Ts, coeffs)

    damping = 0.5;
    alpha0 = 0.24086;

    Torque = Torque_function(theta, activation, damping, alpha0, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6));
    MTF = MTU_function(theta, activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6));
    [LMT, MA] = simplify_pulley_model(coeffs(1), coeffs(2), coeffs(3), theta * 180 / pi);

end
