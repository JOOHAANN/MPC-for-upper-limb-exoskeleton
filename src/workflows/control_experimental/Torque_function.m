function model_torque = Torque_function(activation, damping, alpha0, Ts, L2, L3, R, LoM_nominal, LsT, FoM, L, a, b, c)

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

    MTF(1) = 40; %#play

    % cosa(1)= sqrt(((LM(1)-Ft(1)/k_t)/w)^2/(1+((LM(1)-Ft(1)/k_t)/w)^2));

    alpha(1) = 0.23; %#play

    % dLMT, panation angle, LM, dLM play this!!!

    LM(1) = 0.05;

    LMT(1) = 0.22;

    dLMT(1) = (0.164703 - 0.164693) / 0.007;

    dLM(1) = (0.0376276 - 0.0376376) / 0.007;

    MA(1) = 0.1;

    theta(1) = 1.6;

    dtheta(1) = 0;

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

        ddtheta(k) = MTF(k) * MA(k) / ((a * theta(k)^2 + b * theta(k) + c) * (L^2)) - g / L;

        %% integration and compute Ft

        theta(k + 1) = theta(k) + dtheta(k) * Ts;
        dtheta(k + 1) = dtheta(k) + ddtheta(k) * Ts;
        % LMT(k+1) = sqrt(L3^2 + L2^2 + 2*L2*L3*cos(theta(k+1)));
        [LMT(k + 1), MA(k + 1)] = simplify_pulley_model(L2, L3, R, theta(k + 1));
        dLMT(k + 1) = L2 * L3 * sin(theta(k + 1)) * dtheta(k + 1);

        dFt(k) = k_t(k) * (dLMT(k) + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k))); % handnotes equ (37) and (16)

        MTF(k + 1) = MTF(k) + dFt(k) * Ts;

        % LMT(k+1) = LMT(k) + dLMT(k)*Ts;

        if MTF(k + 1) < 0
            MTF(k + 1) = 0;
        end

        if isnan(MTF(k + 1))
            MTF(k + 1) = MTF(k);
        end

        if MTF(k + 1) > 7000
            MTF(k + 1) = MTF(k);
        end

        t(k + 1) = k * Ts;

    end

    model_torque = MTF(1:end - 1)' .* MA';

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
