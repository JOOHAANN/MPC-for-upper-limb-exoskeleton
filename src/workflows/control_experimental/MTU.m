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

function x_next = MTU(x_current, k, activation, Ts, scale)

    kt_spline = [6.07450471822638e+42    -1.17093401235721e+42   9.78861307802160e+40    -4.68161016668128e+39   1.41111335336212e+38    -2.75023728531764e+36   3.32535388381590e+34    -1.97635568280683e+32   -4.74852712911846e+29   1.68812705982292e+28    -9.30681300424510e+25   -1.69227796153188e+23   3.64048821058250e+21    -7.76789223923884e+18   -4.91004587706162e+16   201470178531589 239829087287.116    -1720265398.11745   -15739.3275509684   6708.90359529906    3.90612013684262];

    damping = 0.5;

    alpha0 = 0.24086;

    PercentageChange = 0.15;

    MTF(k) = x_current(1);
    theta(k) = x_current(2);
    if k > 1
        LT(k - 1) = x_current(3);
    end

    L2 = scale(1);
    L3 = scale(2);
    R = scale(3);
    LoM_nominal = scale(4);
    LsT = scale(5);
    FoM = scale(6);

    w = LoM_nominal * sin(alpha0);

    LMT(k) = sqrt(L3^2 + L2^2 + 2 * L2 * L3 * cos(theta(k)) - R^2);
    %% compute LM for time k

    if k == 1
        k_t(k) = 10 * FoM / LsT;
        dLMT(k) = 0;
    else
        strain(k) = (LT(k - 1) - LsT) / LsT;
        if strain(k) > 0.0199
            strain(k) = 0.0199;
            stiffness(k) = polyval(kt_spline, strain(k));
        else
            stiffness(k) = polyval(kt_spline, strain(k));
        end
        k_t(k) = stiffness(k) * FoM / LsT;
        dLMT(k) = (LMT(k) - LMT(k - 1)) / Ts;
    end

    LoM(k) = LoM_nominal * (1 + PercentageChange * (1 - activation));

    V_max(k) = 10 * LoM(k);

    LT(k) = LsT + MTF(k) / k_t(k);                                                % tendon length
    alpha(k) = acos(sqrt((LMT(k) - LT(k))^2 / (w^2 + (LMT(k) - LT(k))^2)));         % compute the pennation angles for the time k
    LM(k) = (LMT(k) - LT(k)) / cos(alpha(k));                                         % muscle length

    %% compute dLM

    F_L(k) = exp(-17.33 * (abs((LM(k) / LoM(k)) - 1)^3)); % equ (42)
    F_PE(k) = FoM * 3 * 10000 * exp(6 * ((LM(k) / LoM(k)) - 3.2)); % equ (12)
    % F_V(k)=(MTF(k) / cos(alpha(k))-F_PE(k)) / (activation(k) * FoM * exp(-17.33 * abs((LM(k) / LoM(k)) - 1)^3));%solve directly

    %% 1. directive method inverse solving

    %% 2. regression method

    e = [2.2 1];                                                             %  Fitted exactly to the linear part of Fv (k2 in the "Testing_Fv" file
    Z1 = MTF(k) / (FoM * cos(alpha(k))) - F_PE(k) / FoM - activation * F_L(k) * e(2);
    Z2 = activation * F_L(k) * e(1) / (V_max(k)) + damping / V_max(k);
    % dLM(k)=Vw_max * LoM * (0.5+exp(-3*(F_V(k)-0.02))+0.1*exp(12*(F_V(k)-1.5)));

    dLM(k) = Z1 / Z2;

    %% 3. linearized method
    % if k > 1
    % dLM(k)=(LM(k)-LM(k-1))/Ts;
    % end

    dalpha(k) = -dLM(k) / LM(k) * tan(alpha(k));

    % ddLMT(k)= -g/Weight(k)*(Ft(k)+Fac(k) - Weight(k));

    %% integration and compute Ft

    dFt(k) = k_t(k) * (dLMT(k) + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k))); % handnotes equ (37) and (16)

    % dLMT(k+1) = dLMT(k)+ddLMT(k)*Ts;

    % LMT(k+1) = LMT(k) + dLMT(k)*Ts;

    MTF(k + 1) = MTF(k) + dFt(k) * Ts;

    if MTF(k + 1) < 0
        MTF(k + 1) = 0;
    end

    if isnan(MTF(k + 1))
        MTF(k + 1) = MTF(k);
    end

    x_next = [MTF(k + 1), LT(k)];

end
