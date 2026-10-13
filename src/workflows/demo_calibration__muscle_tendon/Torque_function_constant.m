% TORQUE_FUNCTION_CONSTANT Integrate the lumped muscle and exoskeleton dynamics.
% Use prepare_workflow to select matching code and input data.
function X = Torque_function_constant( ...
                                      Tac, theta_real, theta_estimate, activation, activation_front, damping, ...
                                      alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, LMT0)

    % TORQUE_FUNCTION_CONSTANT (clean version without renaming)
    % -------------------------------------------------------------------------
    % Cleaned, structured version of your original function
    % - Keeps all your variable names (LsT, FoM, L1, L2, etc.)
    % - Removes NaN and blow-up risks
    % - Keeps same output format (X with 11 columns)
    % -------------------------------------------------------------------------

    %% ===== Constants =====
    g = 9.81;                 % gravity [m/s^2]
    PercentageChange = 0.15;  % LoM modulation
    vmaxFactor = 10;          % Vmax scaling factor
    initPennation = 0.23;     % initial pennation angle
    initTendonLength = 0.05;  % initial tendon length
    alphaFilterA = 0.92;      % activation smoothing (EMA)
    antiNaNMaxForce = 5e4;    % cap on force
    jointDamping = 10;        % joint viscous damping coefficient
    % zeta = 0.3;               % The ratio of JointDamping

    % PID & control gains
    Kp = 200;
    Kd = 80;
    Ki = 0;
    P = 5;
    D = 3;

    % activation rate limits
    actRateUp = 0.005;
    actRateDown = -0.3;

    % Tendon stiffness polynomial
    kt_poly = [6.07450471822638e+42 -1.17093401235721e+42 9.78861307802160e+40 ...
               -4.68161016668128e+39 1.41111335336212e+38 -2.75023728531764e+36 ...
               3.32535388381590e+34 -1.97635568280683e+32 -4.74852712911846e+29 ...
               1.68812705982292e+28 -9.30681300424510e+25 -1.69227796153188e+23 ...
               3.64048821058250e+21 -7.76789223923884e+18 -4.91004587706162e+16 ...
               201470178531589 239829087287.116 -1720265398.11745 -15739.3275509684 ...
               6708.90359529906 3.90612013684262];

    %% ===== Preallocations =====
    NumberOfData = length(activation);
    theta = zeros(NumberOfData, 1);
    dtheta = zeros(NumberOfData, 1);
    ddtheta = zeros(NumberOfData, 1);
    tau_lim = zeros(NumberOfData, 1);
    % jointDamping = zeros(NumberOfData,1);

    dtheta_real = [0; diff(theta_real) ./ Ts];
    dtheta_estimate = [0; diff(theta_estimate) ./ Ts];
    dLMT = -R * dtheta_real;

    LoM = zeros(NumberOfData, 1);
    alpha = zeros(NumberOfData, 1);
    dalpha = zeros(NumberOfData, 1);
    V_max = zeros(NumberOfData, 1);
    LMT = zeros(NumberOfData, 1);
    LT = zeros(NumberOfData, 1);
    LM = zeros(NumberOfData, 1);
    dLM = zeros(NumberOfData, 1);
    F_L = zeros(NumberOfData, 1);
    F_PE = zeros(NumberOfData, 1);
    T_balance = zeros(NumberOfData, 1);
    F_k_up = zeros(NumberOfData, 1);
    T_front = zeros(NumberOfData, 1);
    strain = zeros(NumberOfData, 1);
    stiffness = zeros(NumberOfData, 1);
    k_t = zeros(NumberOfData, 1);
    MTF = zeros(NumberOfData + 1, 1);
    dFt = zeros(NumberOfData, 1);

    activation = activation(:);
    Tac = Tac(:);

    %% ===== Initial conditions =====
    w = LoM_nominal * sin(alpha0);
    MTF(1) = Torque0 / R;
    alpha(1) = initPennation;
    LT(1) = initTendonLength;
    LMT(1) = LMT0;
    theta(1) = theta0;
    InT = 0;

    %% ===== Main loop =====
    for k = 1:NumberOfData - 1

        % --- tendon stiffness ---
        if k == 1
            k_t(k) = 10 * FoM / LsT;
        else
            strain(k) = (LT(k - 1) - LsT) / LsT;
            if strain(k) > 0.0199
                strain(k) = 0.0199;
            end
            stiffness(k) = polyval(kt_poly, strain(k));
            k_t(k) = stiffness(k) * FoM / LsT;
        end

        % --- muscle-tendon lengths ---
        LMT(k + 1) = LMT(k) + dLMT(k) * Ts;
        LoM(k) = LoM_nominal * (1 + PercentageChange * (1 - activation(k)));
        V_max(k) = vmaxFactor * LoM(k);
        LT(k) = LsT + MTF(k) / k_t(k);
        alpha(k) = acos(sqrt((LMT(k) - LT(k))^2 / (w^2 + (LMT(k) - LT(k))^2)));
        LM(k) = (LMT(k) - LT(k)) / max(cos(alpha(k)), 1e-6);

        % --- force terms ---
        F_L(k) = exp(-17.33 * abs((LM(k) / LoM(k)) - 1)^3);
        F_PE(k) = FoM * 3 * 1e4 * exp(6 * ((LM(k) / LoM(k)) - 3.2));

        % --- dLM ---
        Z1 = MTF(k) / (FoM * max(cos(alpha(k)), 1e-6)) - F_PE(k) / FoM - activation(k) * F_L(k) * 1;
        Z2 = (activation(k) * F_L(k) * 2.2 + damping) / max(V_max(k), 1e-9);
        dLM(k) = Z1 / max(Z2, 1e-9);
        dalpha(k) = -dLM(k) / max(LM(k), 1e-9) * tan(alpha(k));

        % --- balance PID ---
        if activation(k) < 0.05
            err = theta_real(k) - theta(k);
            InT = InT + err * Ts;
            T_balance(k) = Kp * err + Kd * (dtheta_real(k) - dtheta(k)) + Ki * InT;
        else
            % if theta(k) > 1.5
            %     ee = 1 - theta(k);
            %     InT = InT + ee*Ts;
            %     F_k_up(k) = 1*Kp*ee + 1*Kd*(0 - dtheta(k)) + Ki*InT;
            % elseif theta(k) > 1.0
            %     ee = 1 - theta(k);
            %     InT = InT + ee*Ts;
            %     F_k_up(k) = 0.2*Kp*ee + 0.2*Kd*(0 - dtheta(k)) + Ki*InT;
            % end
        end

        % --- activation PD update ---
        if Tac(k) ~= 0
            dactivation = -P * (theta(k) - theta_estimate(k)) - D * (dtheta(k) - dtheta_estimate(k));
            if theta(k) < 1
                dactivation = min(max(dactivation, actRateDown), actRateUp);
            end
            activation(k + 1) = activation(k + 1) + dactivation;
            activation(k + 1) = min(max(activation(k + 1), 0), 1);
        end

        % --- activation smoothing ---
        activation(k + 1) = alphaFilterA * activation(k) + (1 - alphaFilterA) * activation(k + 1);

        % --- load schedule & inertia ---
        [m, I] = local_load_and_inertia(k, L1, L2, M);

        % % --- JointDamping ---
        %
        % jointDamping(k) = 2 * zeta * sqrt( I * k_t(k) * R ^ 2 );

        % --- Soft limits: tau_lim ~ 0 inside, active only outside ---
        theta_min = deg2rad(0);
        theta_max = deg2rad(95);

        a  = 30;     % 开关“硬度”（越大越接近硬开关：20–50）
        k_wall = 1200;    % 软墙刚度（线性区斜率），按需要调大/小
        Fmax   = 900;    % 软墙力上限（防爆力；≥系统可能的外推峰值）
        c_wall = 20;     % 墙面阻尼（只在越界且向外运动时生效）

        % 平滑Heaviside（开关），区间内≈0，越界≈1
        H  = @(x) 0.5 * (1 + tanh(a * x));

        % 越界激活（上/下侧）
        actUp   = H(theta(k) - theta_max);   % q>q_max 时≈1，否则≈0
        actDown = H(theta_min - theta(k));   % q<q_min 时≈1，否则≈0

        % 渗透量：区间内=0；越界=穿透深度（用 ReLU 更“干净”）
        penU = max(theta(k) - theta_max, 0);   % 上侧穿透
        penL = max(theta_min - theta(k), 0);   % 下侧穿透

        % 线性回拉 + 平滑饱和，仅在对应侧激活
        tau_up   = -Fmax * tanh((k_wall / Fmax) * penU) * actUp;
        tau_down = +Fmax * tanh((k_wall / Fmax) * penL) * actDown;

        % 单向墙面阻尼：只阻尼“继续往外推”的速度
        v_out_up   = max(dtheta(k),  0);   % 朝上越界方向的速度分量
        v_out_down = max(-dtheta(k),  0);   % 朝下越界方向的速度分量
        tau_damp   = -c_wall * (v_out_up * actUp - v_out_down * actDown);

        % 合力（若你还想保留原来的通用粘性项，可另外再加）
        tau_lim(k) = tau_up + tau_down + tau_damp;

        % --- dynamics ---
        tauGrav = -m * g * L2 * cos(theta(k) - phi) - M * g * L1 * cos(theta(k));
        ddtheta(k) = (-jointDamping * dtheta(k) + F_k_up(k) + MTF(k) * R + Tac(k) + ...
                      T_balance(k) + tauGrav + tau_lim(k)) / max(I, 1e-9);
        dtheta(k + 1) = dtheta(k) + ddtheta(k) * Ts;
        theta(k + 1) = theta(k) + dtheta(k) * Ts;

        % --- muscle-tendon force update ---
        dFt(k) = k_t(k) * (dLMT(k) + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k)));
        MTF(k + 1) = MTF(k) + dFt(k) * Ts;

        % --- safety ---
        if isnan(MTF(k + 1))
            MTF(k + 1) = MTF(k);
        end
        if MTF(k + 1) < 0
            MTF(k + 1) = 0;
        end
        if MTF(k + 1) > antiNaNMaxForce
            MTF(k + 1) = antiNaNMaxForce;
        end
    end

    %% ===== Output =====
    model_torque = MTF(1:end - 1) * R;
    X = zeros(NumberOfData, 11);
    X(:, 1) = model_torque;
    X(:, 2) = theta;
    X(:, 3) = F_L;
    X(:, 4) = LMT;
    X(:, 5) = k_t;
    X(:, 6) = tau_lim;
    X(:, 7) = activation;
    X(:, 8) = F_PE;
    X(:, 9) = Tac;
    X(:, 10) = T_balance;
    X(:, 11) = MTF(1:end - 1) - F_PE;

end

%% ===== Helper =====
function [mLoad, I] = local_load_and_inertia(k, L1, L2, M)
    if      k > 9050  && k < 12537
        m = 5;
    elseif  k > 22662 && k < 26159
        m = 10;
    elseif  k > 35678 && k < 37606
        m = 15;
    else
        m = 0;
    end
    mLoad = m;
    I = m * L2^2 + (4 * M * L1^2) / 3;
end
