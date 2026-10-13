function x = Torque_function_constant(m, Tac, theta_real, theta_estimate, activation, activation_front, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, LMT0, front_a, front_b, front_c)
    % TORQUE_FUNCTION_CONSTANT Integrate the lumped muscle and exoskeleton dynamics.
    % Use prepare_workflow to select matching code and input data.

    PercentageChange = 0.15;

    g = 9.81;

    NumberOfData = length(activation);

    %% Preallocations
    theta = zeros(NumberOfData, 1)';
    dtheta = zeros(NumberOfData, 1)';
    ddtheta = zeros(NumberOfData, 1)';
    dtheta_real = zeros(NumberOfData, 1)';
    dtheta_estimate = zeros(NumberOfData, 1)';

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
    T_balance = zeros(NumberOfData, 1)';
    T_front = zeros(NumberOfData, 1)';

    MTF = zeros(NumberOfData + 1, 1)';
    dFt = zeros(NumberOfData, 1)';
    % Tac = zeros(NumberOfData,1)';
    strain = zeros(NumberOfData, 1)';

    kt_spline = [6.07450471822638e+42   -1.17093401235721e+42   9.78861307802160e+40    -4.68161016668128e+39   1.41111335336212e+38    -2.75023728531764e+36   3.32535388381590e+34    -1.97635568280683e+32   -4.74852712911846e+29   1.68812705982292e+28    -9.30681300424510e+25   -1.69227796153188e+23   3.64048821058250e+21    -7.76789223923884e+18   -4.91004587706162e+16   201470178531589 239829087287.116    -1720265398.11745   -15739.3275509684   6708.90359529906    3.90612013684262];

    w = LoM_nominal * sin(alpha0);

    %% compute the pennation angle for the first time

    MTF(1) = Torque0 / R; %#play

    % cosa(1)= sqrt(((LM(1)-Ft(1)/k_t)/w)^2/(1+((LM(1)-Ft(1)/k_t)/w)^2));

    alpha(1) = 0.23; %#play

    LT(1) = 0.05;

    LMT(1) = LMT0;

    dLMT(1) = 0;

    dLM(1) = 0;

    MA(1) = R;

    theta(1) = theta0;

    dtheta(1) = 0;

    a = 0.92;

    c = 10;

    Kp = 200;

    Kd = 80;

    Ki = 0;

    InT = 0;

    % P = 0.3;
    %
    % D = 0.05;

    P = 5;

    D = 3;

    dtheta_real = [0; diff(theta_real) ./ Ts];

    dtheta_estimate = [0; diff(theta_estimate) ./ Ts];

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

        %% 2. linearized method

        e = [2.2 1];                                                             %  Fitted exactly to the linear part of Fv (k2 in the "Testing_Fv" file
        Z1 = MTF(k) / (FoM * cos(alpha(k))) - F_PE(k) / FoM - activation(k) * F_L(k) * e(2);
        Z2 = activation(k) * F_L(k) * e(1) / (V_max(k)) + damping / V_max(k);
        % dLM(k)=Vw_max * LoM * (0.5+exp(-3*(F_V(k)-0.02))+0.1*exp(12*(F_V(k)-1.5)));
        dLM(k) = Z1 / Z2;

        dalpha(k) = -dLM(k) / LM(k) * tan(alpha(k));

        % ddtheta(k) = 0;
        MA(k) = R;

        % if theta(k) > 1.7
        %     % activation(k+1) = 0;
        %     k_lock_up = 2000;
        % else
        %     k_lock_up = 0;
        % end
        % if (theta(k) > 1.2)  && (activation(k) > 0.05)
        %     ee = 1.2 - theta(k);
        %     InT = InT + ee*Ts;
        %     F_k_up(k) = Kp*(ee) + Kd*(0 - dtheta(k)) + Ki*InT;
        % end

        %
        % F_k_down(k) =  k_lock_down*(abs(theta(k)-(-0.5))) + k_lock_down*(abs(dtheta(k) - 0));
        % F_k_up(k) =  k_lock_up*(abs(theta(k) - 1.7)) + 0*k_lock_up*(abs(dtheta(k) - 0));
        % F_k(k) =  0;

        % %%%%%%%%%%%%%%%%%%%PID using to track the real theta%%%%%%%%%%%%%%%%%%%%%%%

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

        % %%%%%%%%%%%%%%%%%%%PID using to track the real theta%%%%%%%%%%%%%%%%%%%%%%%

        % %%%%%%%%%%%%%%%%%%PID using to track the real theta%%%%%%%%%%%%%%%%%%%%%%%
        if Tac(k) ~= 0
            dactivation = -P * (theta(k) - theta_estimate(k)) - D * (dtheta(k) - dtheta_estimate(k));
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
        % if activation(k+1) - activation(k) > 0.01
        %     activation(k+1) = activation(k) + 0.01;
        % end
        %% the activation filter to smooth the data

        % %%%%%%%%%%%%%%%%%%PID using to track the real theta%%%%%%%%%%%%%%%%%%%%%%%

        % % %%%%%%%step4%%%%%%%%Place the loads on the model%%%%%%%%%%%%%%%%%%%%%%%

        if k > 9050 && k < 12537
            m = 5;
            I = m * L2^2 + (4 * M * L1^2) / 3;
        elseif k > 22662 && k < 26159
            m = 10;
            I = m * L2^2 + (4 * M * L1^2) / 3;
        elseif k > 37966 && k < 40544
            m = 15;
            I = m * L2^2 + (4 * M * L1^2) / 3;
        else
            m = 0;
            I = m * L2^2 + (4 * M * L1^2) / 3;
        end

        % % %%%%%%%step4%%%%%%%%Place the loads on the model%%%%%%%%%%%%%%%%%%%%%%%

        % % %%%%%%%%%%%%%%%%%%Front rect muscle torque%%%%%%%%%%%%%%%%%%%%%%%

        T_front(k) = front_a * (theta(k) - 0.2)^2 + front_b * (theta(k) - 0.2) + front_c;

        % % %%%%%%%%%%%%%%%%%%Front rect muscle torque%%%%%%%%%%%%%%%%%%%%%%%

        % if  Torque_real_1(k) > 80
        %     ddtheta(k) = (-c*dtheta(k) - F_k_up(k) + F_k_down(k) + MTF(k)*MA(k) + Tac(k) + T_balance(k) - T_front(k) - m*g*L2*cos(theta(k)-phi)-M*g*L1*cos(theta(k))) / I;
        % % ddtheta(k) = (MTF(k)*MA(k) + Tac - (M+m)*g*L1*theta(k)^2 + L2*theta(k)) / phi;
        % else
        %     ddtheta(k) = (-c*dtheta(k) - F_k_up(k) + F_k_down(k) + MTF(k)*MA(k) + 0 + T_balance(k) - T_front(k) - m*g*L2*cos(theta(k)-phi)-M*g*L1*cos(theta(k))) / I;
        % end
        ddtheta(k) = (-c * dtheta(k) + F_k_up(k) + F_k_down(k) + MTF(k) * MA(k) + Tac(k) + T_balance(k) - T_front(k) - m * g * L2 * cos(theta(k) - phi) - M * g * L1 * cos(theta(k))) / I;

        dtheta(k + 1) = dtheta(k) + ddtheta(k) * Ts;
        theta(k + 1) = theta(k) + dtheta(k) * Ts;

        %% integration and compute Ft

        LMT(k + 1) = LMT(k) + dLMT(k) * Ts;

        dLMT(k + 1) = -R * dtheta(k + 1);

        dFt(k) = k_t(k) * (dLMT(k) + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k))); % handnotes equ (37) and (16)

        MTF(k + 1) = MTF(k) + dFt(k) * Ts;

        % LMT(k+1) = LMT(k) + dLMT(k)*Ts;

        if MTF(k + 1) < -0
            MTF(k + 1) = -0;
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
    x(:, 3) = F_L;
    x(:, 4) = LMT;
    x(:, 5) = F_k_up;
    x(:, 6) = T_front;
    x(:, 7) = activation;
    x(:, 8) = F_PE;
    x(:, 9) = Tac;
    x(:, 10) = T_balance;
    x(:, 11) = MTF(1:end - 1) - F_PE;
