function x = Torque_function_constant(m, Tac, Real_theta, activation, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, LMT0)
    % TORQUE_FUNCTION_CONSTANT Integrate the lumped muscle and exoskeleton dynamics.
    % Use prepare_workflow to select matching code and input data.

    PercentageChange = 0.15;

    g = 9.81;

    NumberOfData = length(activation);

    %% Preallocations
    theta = zeros(NumberOfData, 1)';
    dtheta = zeros(NumberOfData, 1)';
    ddtheta = zeros(NumberOfData, 1)';

    LoM = zeros(NumberOfData, 1)';
    alpha = zeros(NumberOfData, 1)';
    dalpha = zeros(NumberOfData, 1)';
    V_max = zeros(NumberOfData, 1)';
    LMT = zeros(NumberOfData, 1)';
    MA = zeros(NumberOfData, 1)';
    LM = zeros(NumberOfData, 1)';
    dLM = zeros(NumberOfData, 1)';
    LT = zeros(NumberOfData, 1)';
    dLMT = zeros(NumberOfData, 1)';
    k_t = zeros(NumberOfData, 1)';
    stiffness = zeros(NumberOfData, 1)';
    F_L = zeros(NumberOfData, 1)';
    F_PE = zeros(NumberOfData, 1)';

    t = zeros(NumberOfData, 1)';

    MTF = zeros(NumberOfData + 1, 1)';
    dFt = zeros(NumberOfData, 1)';
    % Tac = zeros(NumberOfData,1)';
    strain = zeros(NumberOfData, 1)';

    kt_spline = [6.07450471822638e+42   -1.17093401235721e+42   9.78861307802160e+40    -4.68161016668128e+39   1.41111335336212e+38    -2.75023728531764e+36   3.32535388381590e+34    -1.97635568280683e+32   -4.74852712911846e+29   1.68812705982292e+28    -9.30681300424510e+25   -1.69227796153188e+23   3.64048821058250e+21    -7.76789223923884e+18   -4.91004587706162e+16   201470178531589 239829087287.116    -1720265398.11745   -15739.3275509684   6708.90359529906    3.90612013684262];

    w = LoM_nominal * sin(alpha0);

    % %% compute the LMT of the pulley model
    %
    %
    % LMT = sqrt(L3^2 + L2^2 + 2*L2*L3*cos(theta)-R^2);
    %
    % AC = L3^2 + L2^2 + 2*L2*L3*cos(theta);
    %
    % C1 = (AC + L2^2 - L3^2) ./ (2 * L2 * sqrt(AC));
    %
    % Sa1a2 = sqrt(1-C1.^2) .* LMT./sqrt(AC) + C1 .* R./sqrt(AC);
    %
    % MA = Sa1a2 * L2;

    %% compute the pennation angle for the first time

    MTF(1) = 100; %#play

    % cosa(1)= sqrt(((LM(1)-Ft(1)/k_t)/w)^2/(1+((LM(1)-Ft(1)/k_t)/w)^2));

    alpha(1) = 0.23; %#play

    % dLMT, panation angle, LM, dLM play this!!!
    LT(1) = 0.05;

    LMT(1) = LMT0;

    dLMT(1) = 0;

    dLM(1) = 0;

    MA(1) = R;

    theta(1) = 0.66;

    dtheta(1) = 0;

    L = 0.3;

    a = 0.99;

    I = m * L2^2 + (4 * M * L1^2) / 3;

    for k = 1:(NumberOfData - 1)

        %% compute LM for time k

        if k == 1
            k_t(k) = 10 * FoM / LsT;
        else
            strain(k) = (LT(k - 1) - LsT) / LsT;
            if strain(k) > 0.0199
                strain(k) = 0.0199;
                stiffness(k) = polyval(kt_spline, strain(k));
            else
                stiffness(k) = polyval(kt_spline, strain(k));
            end
            k_t(k) = stiffness(k) * FoM / LsT;
        end

        if k > 2500
            activation(k) = 0;
        end
        LoM(k) = LoM_nominal * (1 + PercentageChange * (1 - activation(k)));

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
        Z1 = MTF(k) / (FoM * cos(alpha(k))) - F_PE(k) / FoM - activation(k) * F_L(k) * e(2);
        Z2 = activation(k) * F_L(k) * e(1) / (V_max(k)) + damping / V_max(k);
        % dLM(k)=Vw_max * LoM * (0.5+exp(-3*(F_V(k)-0.02))+0.1*exp(12*(F_V(k)-1.5)));

        dLM(k) = Z1 / Z2;

        %% 3. linearized method
        % if k > 1
        % dLM(k)=(LM(k)-LM(k-1))/Ts;
        % end

        dalpha(k) = -dLM(k) / LM(k) * tan(alpha(k));

        % ddtheta(k) = 0;
        MA(k) = R;
        % if MTF(k)*MA(k) > 50
        %     ddtheta(k) = (MTF(k)*MA(k) + 50 - M) / I;
        % else
        %     ddtheta(k) = (MTF(k)*MA(k) + 0 - M) / I;
        % end
        ddtheta(k) = (MTF(k) * MA(k) + Tac - m * g * L2 * cos(theta(k) - phi) - M * g * L1 * cos(theta(k))) / I;
        dtheta(k + 1) = dtheta(k) + ddtheta(k) * Ts;
        theta(k + 1) = theta(k) + dtheta(k) * Ts;

        %% integration and compute Ft

        LMT(k + 1) = LMT(k) + dLMT(k) * Ts;

        dLMT(k + 1) = -R * dtheta(k + 1);

        dFt(k) = k_t(k) * (dLMT(k) + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k))); % handnotes equ (37) and (16)

        MTF(k + 1) = MTF(k) + dFt(k) * Ts;

        MTF(k + 1) =  a * MTF(k + 1) + (1 - a) * MTF(k);

        % LMT(k+1) = LMT(k) + dLMT(k)*Ts;

        if MTF(k + 1) < -5000
            MTF(k + 1) = -5000;
        end

        if isnan(MTF(k + 1))
            MTF(k + 1) = MTF(k);
        end

        if MTF(k + 1) > 50000
            MTF(k + 1) = MTF(k);
        end

        t(k + 1) = k * Ts;

    end

    model_torque = MTF(1:end - 1)' .* MA';

    LM_LoM = LM ./ mean(LoM);

    x(:, 1) = model_torque;
    x(:, 2) = theta;
    x(:, 3) = LMT;
    x(:, 4) = ddtheta;
    x(:, 5) = LM;
    x(:, 6) = dLMT;
    x(:, 7) = LoM;
    x(:, 8) = F_PE;
    x(:, 9) = MTF(1:end - 1);
    x(:, 10) = LM_LoM;

    % %% plot figure
    %
    % MTF = MTF';
    %
    % figure;
    % plot(t(1:750),MTF(1:750));
    % ylim([-40, 90]);
    % hold on
    % plot(t(1:750), MTF_real(1:750));
    % plot(t(1:750), activation(1:750));
    % xlabel('t');
    % ylabel('Force');
    % ylim([-200, 1000]);
    % legend('MTF', 'F_desired','activation');
    % title('Ft(t) , F desired force, and 1000*activation');
    % grid on;
