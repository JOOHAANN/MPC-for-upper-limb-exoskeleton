% % 定义非线性系统
% function x_next = nonlinear_system(x, u, Ts)
%     % 非线性系统的动力学方程，离散化
%     MTF = x(1);
%     theta = x(2);
%     angle_velocity = x(3);
%     theta = x(4);
%     x_next = [x1 + Ts * x2;
%               x2 + Ts * u * (1 - x2^2)];
% end

function [MTF_next, theta_next, dtheta_next] = nonlinear_system(k_t, InitialConditions, activation, u_Tac, Ts, coeffs)

    L2 = coeffs(1);
    L3 = coeffs(2);
    R = coeffs(3);
    LoM_nominal = coeffs(4);
    LsT = coeffs(5);
    FoM = coeffs(6);
    I = coeffs(7);
    G = coeffs(8);

    PercentageChange = 0.15;

    g = 9.81;

    kt_spline = [6.07450471822638e+42   -1.17093401235721e+42   9.78861307802160e+40    -4.68161016668128e+39   1.41111335336212e+38    -2.75023728531764e+36   3.32535388381590e+34    -1.97635568280683e+32   -4.74852712911846e+29   1.68812705982292e+28    -9.30681300424510e+25   -1.69227796153188e+23   3.64048821058250e+21    -7.76789223923884e+18   -4.91004587706162e+16   201470178531589 239829087287.116    -1720265398.11745   -15739.3275509684   6708.90359529906    3.90612013684262];

    damping = 0.5;

    alpha0 = 0.24086;

    w = LoM_nominal * sin(alpha0);

    %% compute the pennation angle for the first time

    MTF(1) = InitialConditions(1); % state x1

    theta(1) = InitialConditions(2); % state x2

    dtheta(1) = InitialConditions(3); % state x3

    [LMT(1), MA(1)] = simplify_pulley_model(L2, L3, R, theta(1));

    %% compute LM for time k

    % strain(1)= (LT(1) - LsT)/LsT;
    % if strain(1) > 0.0199
    %     strain(1) = 0.0199;
    %     stiffness(1)= polyval(kt_spline,strain(1));
    % else
    %     stiffness(1)= polyval(kt_spline,strain(1));
    % end
    % k_t(1)= stiffness(k)*FoM/LsT;
    k_t(1) = k_t;

    LoM(1) = LoM_nominal * (1 + PercentageChange * (1 - activation(1)));

    V_max(1) = 10 * LoM(1);

    LT(1) = LsT + MTF(1) / k_t(1);                                                % tendon length
    alpha(1) = acos(sqrt((LMT(1) - LT(1))^2 / (w^2 + (LMT(1) - LT(1))^2)));      % compute the pennation angles for the time k
    LM(1) = (LMT(1) - LT(1)) / cos(alpha(1));                                   % muscle length

    %% compute dLM

    F_L(1) = exp(-17.33 * (abs((LM(1) / LoM(1)) - 1)^3)); % equ (42)
    F_PE(1) = FoM * 3 * 10000 * exp(6 * ((LM(1) / LoM(1)) - 3.2)); % equ (12)

    e = [2.2 1];                                                                 %  Fitted exactly to the linear part of Fv (k2 in the "Testing_Fv" file
    Z1 = MTF(1) / (FoM * cos(alpha(1))) - F_PE(1) / FoM - activation(1) * F_L(1) * e(2);
    Z2 = activation(1) * F_L(1) * e(1) / (V_max(1)) + damping / V_max(1);

    dLM(1) = Z1 / Z2;

    dalpha(1) = -dLM(1) / LM(1) * tan(alpha(1));

    ddtheta(1) = (MTF(1) * MA(1) + u_Tac(1) - G) / I;

    %% integration and compute Ft

    theta(2) = theta(1) + dtheta(1) * Ts;
    dtheta(2) = dtheta(1) + ddtheta(1) * Ts;
    % LMT(k+1) = sqrt(L3^2 + L2^2 + 2*L2*L3*cos(theta(k+1)));
    % [LMT(2), MA(2)] = simplify_pulley_model(L2, L3, R, theta(2));              %the state of LMT and MA

    dLMT(2) = L2 * L3 * sin(theta(2)) * dtheta(2);

    dFt(1) = k_t(1) * (dLMT(1) + LM(1) * sin(alpha(1)) * dalpha(1) - dLM(1) * cos(alpha(1))); % handnotes equ (37) and (16)

    MTF(2) = MTF(1) + dFt(1) * Ts;

    % LMT(k+1) = LMT(k) + dLMT(k)*Ts;

    if MTF(2) < 0
        MTF(2) = 0;
    end

    if isnan(MTF(2))
        MTF(2) = MTF(1);
    end

    if MTF(2) > 7000
        MTF(2) = MTF(1);
    end

    MTF_next = MTF(2);
    theta_next = theta(2);
    dtheta_next = dtheta(2);

end
