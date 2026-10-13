function x = nonlinear_system_period(m, Tac, InitialConditions, activation, theta_real, dtheta_real, damping, alpha0, Ts, coeffs)

    PercentageChange = 0.15;

    g = 9.81;

    NumberOfData = length(activation);

    R = coeffs(1);
    LoM_nominal = coeffs(2);
    LsT = coeffs(3);
    FoM = coeffs(4);
    M = coeffs(5);
    L1 = coeffs(6);
    L2 = coeffs(7);
    phi = coeffs(8);
    front_a = coeffs(12);
    front_b = coeffs(13);
    front_c = coeffs(14);
    c = 10;
    Kp = 200;
    Kd = 80;
    Ki = 0;
    InT = 0;
    P = 1;
    D = 0.6;
    I = m * L2^2 + (4 * M * L1^2) / 3;
    a = 0.9;

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
    F_k_up = zeros(NumberOfData, 1)';
    F_k_down = zeros(NumberOfData, 1)';
    T_front = zeros(NumberOfData, 1)';
    T_balance = zeros(NumberOfData, 1)';

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

    MTF(1) = InitialConditions(1); % state x1
    theta(1) = InitialConditions(2); % state x2
    dtheta(1) = InitialConditions(3); % state x3
    LMT(1) = InitialConditions(4); % state x4
    dLMT(1) = InitialConditions(5); % state x5
    LT(1) = InitialConditions(6); % state x6

    for k = 1:(NumberOfData - 1)

        %% compute LM for time k

        strain(k) = (LT(k) - LsT) / LsT;
        if strain(k) > 0.0199
            strain(k) = 0.0199;
            stiffness(k) = polyval(kt_spline, strain(k));
        else
            stiffness(k) = polyval(kt_spline, strain(k));
        end
        k_t(k) = stiffness(k) * FoM / LsT;

        LoM(k) = LoM_nominal * (1 + PercentageChange * (1 - activation(k)));

        V_max(k) = 10 * LoM(k);

        LT(k + 1) = LsT + MTF(k) / k_t(k);                                                % tendon length
        alpha(k) = acos(sqrt((LMT(k) - LT(k + 1))^2 / (w^2 + (LMT(k) - LT(k + 1))^2)));         % compute the pennation angles for the time k
        LM(k) = (LMT(k) - LT(k + 1)) / cos(alpha(k));                                         % muscle length

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

        % if theta(k) > 1.7
        %     % activation(k+1) = 0;
        %     k_lock_up = 2000;
        % else
        %     k_lock_up = 0;
        % end
        %
        % F_k_up(k) =  k_lock_up*(abs(theta(k) - 1.7)) + 0.5*k_lock_up*(abs(dtheta(k) - 0));
        % F_k(k) =  0;

        % %%%%%%%%%%%%%%%%%%%PID using to track the theta = 1.5%%%%%%%%%%%%%%%%%%%%%%%

        % if (k < 7757) || (k > 13853 && k < 20956) || (k > 27639 && k < 36908) || (k > 41734) || (theta(k) > 1.5)
        % if k < 7757 || k > 13853
        % || (theta(k) > 1.5)
        if activation(k) < 0.05
            err = theta_real(k) - theta(k);
            InT = InT + err * Ts;
            T_balance(k) = Kp * (err) + Kd * (dtheta_real(k) - dtheta(k)) + Ki * InT;
        else
            if (theta(k) > 1)  && (activation(k) > 0.05)
                ee = 1 - theta(k);
                InT = InT + ee * Ts;
                F_k_up(k) = 0.2 * Kp * (ee) + 0.2 * Kd * (0 - dtheta(k)) + Ki * InT;
            elseif (theta(k) > 1.5)  && (activation(k) > 0.05)
                ee = 1 - theta(k);
                InT = InT + ee * Ts;
                F_k_up(k) = 1 * Kp * (ee) + 1 * Kd * (0 - dtheta(k)) + Ki * InT;
            end
        end

        % %%%%%%%%%%%%%%%%%%%PID using to track the theta = 1.5%%%%%%%%%%%%%%%%%%%%%%%

        % %%%%%%%%%%%%%%%%%%PID using to track the real theta when Tac(k) ~= 0%%%%%%%%%%%%%%%%%%%%%%%
        if Tac(k) ~= 0
            dactivation = -P * (theta(k) - theta_real(k)) - D * (dtheta(k) - dtheta_real(k));
            if theta(k) < 1
                if dactivation > 0.005
                    dactivation = 0.005;
                elseif dactivation < -0.3
                    dactivation = -0.3;
                end
                activation(k + 1) = activation(k + 1) + dactivation;
            else
                activation(k + 1) = activation(k + 1) + dactivation;
            end
            % if (activation(k+1) - activation(k))/Ts > 50
            %     activation(k+1) = 50*Ts + activation(k);
            % end

            % elseif (activation(k+1) - activation(k))/Ts < -0.5
            %     activation(k+1) = -0.5*Ts + activation(k);
            if activation(k + 1) < 0
                activation(k + 1) = 0;
            elseif activation(k + 1) > 1
                activation(k + 1) = 1;
            end
        end

        %% the activation filter to smooth the data
        activation(k + 1) = a * activation(k) + (1 - a) * activation(k + 1);
        %% the activation filter to smooth the data

        % %%%%%%%%%%%%%%%%%%PID using to track the real theta when Tac(k) ~= 0%%%%%%%%%%%%%%%%%%%%%%%
        % % %%%%%%%%%%%%%%%%%%Front rect muscle torque%%%%%%%%%%%%%%%%%%%%%%%

        T_front(k) = front_a * (theta(k) - 0.5)^2 + front_b * (theta(k) - 0.5) + front_c;

        % % %%%%%%%%%%%%%%%%%%Front rect muscle torque%%%%%%%%%%%%%%%%%%%%%%%

        ddtheta(k) = (-c * dtheta(k) + F_k_up(k) + F_k_down(k) + MTF(k) * MA(k) + Tac(k) + T_balance(k) - T_front(k) - m * g * L2 * cos(theta(k) - phi) - M * g * L1 * cos(theta(k))) / I;
        % ddtheta(k) = (MTF(k)*MA(k) + Tac - (M+m)*g*L1*theta(k)^2 + L2*theta(k)) / phi;

        dtheta(k + 1) = dtheta(k) + ddtheta(k) * Ts;
        theta(k + 1) = theta(k) + dtheta(k) * Ts;

        %% integration and compute Ft

        LMT(k + 1) = LMT(k) + dLMT(k) * Ts;

        dLMT(k + 1) = -R * dtheta(k + 1);

        dFt(k) = k_t(k) * (dLMT(k) + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k))); % handnotes equ (37) and (16)

        MTF(k + 1) = MTF(k) + dFt(k) * Ts;

        if MTF(k + 1) < 0
            MTF(k + 1) = 0;
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

    x(:, 1) = model_torque(1:end - 1);
    x(:, 2) = dtheta(1:end - 1);
    x(:, 3) = theta(1:end - 1);

end

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
