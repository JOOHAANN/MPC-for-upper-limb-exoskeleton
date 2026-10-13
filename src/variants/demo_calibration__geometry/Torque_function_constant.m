function x = Torque_function_constant(Tac, theta_real, activation, Ts, M, L1, L2, phi, Torque0, LMT0)
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
    model_torque = zeros(NumberOfData, 1)';

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

    %% compute the pennation angle for the first time

    alpha(1) = 0.23; %#play

    LT(1) = 0.05;

    LMT(1) = LMT0;

    dLMT(1) = 0;

    dLM(1) = 0;

    dtheta(1) = 0;

    a = 0.92;

    c = 10;

    dtheta_real = [0; diff(theta_real) ./ Ts];

    ddtheta_real = [0; diff(dtheta_real) ./ Ts];

    for k = 1:(NumberOfData - 1)

        %% compute LM for time k
        % % %%%%%%%step4%%%%%%%%Place the loads on the model%%%%%%%%%%%%%%%%%%%%%%%

        if k > 9050 && k < 12537
            m = 5;
            I = m * L2^2 + (4 * M * L1^2) / 3;
        elseif k > 22662 && k < 26159
            m = 10;
            I = m * L2^2 + (4 * M * L1^2) / 3;
        elseif k > 35677 && k < 37600
            m = 15;
            I = m * L2^2 + (4 * M * L1^2) / 3;
        else
            m = 0;
            I = m * L2^2 + (4 * M * L1^2) / 3;
        end

        % if k > 2000 && k < 17000
        %     m = 5;
        %     I = m*L2^2 + (4*M*L1^2)/3;
        % elseif k > 17000 && k < 32000
        %     m = 10;
        %     I = m*L2^2 + (4*M*L1^2)/3;
        % elseif k > 32000 && k < 45000
        %     m = 15;
        %     I = m*L2^2 + (4*M*L1^2)/3;
        % else
        %     m = 0;
        %     I = m*L2^2 + (4*M*L1^2)/3;
        % end

        % % %%%%%%%step4%%%%%%%%Place the loads on the model%%%%%%%%%%%%%%%%%%%%%%%

        model_torque(k) = ddtheta_real(k) * I + c * dtheta_real(k) - Tac(k) + m * g * L2 * cos(theta_real(k) - phi) + M * g * L1 * cos(theta_real(k));

        %% integration and compute Ft

        % LMT(k+1) = LMT(k) + dLMT(k)*Ts;

    end

    LM_LoM = LM ./ mean(LoM);

    x(:, 1) = model_torque;
    x(:, 2) = theta_real;
    x(:, 3) = dtheta_real;
    x(:, 4) = ddtheta_real;
    x(:, 5) = F_k_up;
    x(:, 6) = T_front;
    x(:, 7) = activation;
    x(:, 8) = F_PE;
    x(:, 9) = Tac;
    x(:, 10) = T_balance;
    x(:, 11) = MTF(1:end - 1) - F_PE;
