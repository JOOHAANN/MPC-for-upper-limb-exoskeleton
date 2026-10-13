% TORQUE_FUNCTION_CONSTANT Integrate the lumped muscle and exoskeleton dynamics.
% Use prepare_workflow to select matching code and input data.
function X = Torque_function_constant( ...
                                      Tac, theta_real, theta_estimate, activation, activation_front, damping, ...
                                      alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, zeta, LMT0)

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
    % jointDamping = 10;        % joint viscous damping coefficient
    % zeta = 0.1;               % The ratio of JointDamping

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
    jointDamping = zeros(NumberOfData, 1);

    dtheta_real = [0; diff(theta_real) ./ Ts];
    dtheta_estimate = [0; diff(theta_estimate) ./ Ts];
    % dLMT = -R * dtheta_real;

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
    dLMT = zeros(NumberOfData, 1);

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
            strain(k) = min(max(strain(k), 0), 0.0199);
            stiffness(k) = max(polyval(kt_poly, strain(k)), 0);
            k_t(k) = max(stiffness(k) * FoM / LsT, 1e-6);
        end

        % --- muscle-tendon lengths ---
        % LMT(k+1) = LMT(k) + dLMT(k) * Ts;
        LoM(k) = LoM_nominal * (1 + PercentageChange * (1 - activation(k)));
        V_max(k) = vmaxFactor * LoM(k);
        LT(k) = LsT + MTF(k) / k_t(k);
        tmp = (LMT(k) - LT(k));
        arg = (tmp^2) / (w^2 + tmp^2 + 1e-12);
        arg = min(max(arg, 0), 1);     % 关键：钳制到 [0,1]
        alpha(k) = acos(sqrt(arg));

        LM(k) = tmp / max(cos(alpha(k)), 1e-6);

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

        % --- JointDamping ---

        jointDamping(k) = 2 * zeta * sqrt(I * k_t(k) * R^2);

        % ====== (NEW) inner iteration for implicit-ish update ======
        nIt = 3;  % 2~4 generally enough
        theta_next = theta(k);
        dtheta_next = dtheta(k);
        LMT_next   = LMT(k);
        MTF_next   = MTF(k);

        for it = 1:nIt

            % --- kinematics: use predicted states in this iteration ---
            % activation -> LoM, Vmax (you can keep using activation(k) here)
            LoM(k)   = LoM_nominal * (1 + PercentageChange * (1 - activation(k)));
            V_max(k) = vmaxFactor * LoM(k);

            % use predicted LMT/MTF to compute tendon length etc.
            % tendon length from force & stiffness at step k (k_t(k) computed already)
            LT(k) = LsT + MTF_next / max(k_t(k), 1e-9);

            % pennation and fiber length from predicted LMT_next and LT(k)
            % (protect against negative/zero geometry)
            tmp = (LMT_next - LT(k));
            arg = (tmp^2) / (w^2 + tmp^2 + 1e-12);
            arg = min(max(arg, 0), 1);     % 关键：防 acos 复数
            alpha(k) = acos(sqrt(arg));

            cA = max(cos(alpha(k)), 1e-6);
            LM(k) = tmp / cA;

            % --- force terms (same as your code) ---
            F_L(k)  = exp(-17.33 * abs((LM(k) / max(LoM(k), 1e-9)) - 1)^3);
            F_PE(k) = FoM * 3 * 1e4 * exp(6 * ((LM(k) / max(LoM(k), 1e-9)) - 3.2));

            % --- dLM, dalpha (same structure, but based on predicted MTF_next) ---
            Z1 = MTF_next / (FoM * cA) - F_PE(k) / FoM - activation(k) * F_L(k) * 1;
            Z2 = (activation(k) * F_L(k) * 2.2 + damping) / max(V_max(k), 1e-9);
            dLM(k) = Z1 / max(Z2, 1e-9);
            dalpha(k) = -dLM(k) / max(LM(k), 1e-9) * tan(alpha(k));

            % --- dynamics: implicit damping + semi-implicit theta update ---
            % compute tau_lim using predicted theta_next/dtheta_next
            tau_lim_next = local_tau_lim(theta_next, dtheta_next);

            % gravity torque using predicted theta_next
            tauGrav_next = -m * g * L2 * cos(theta_next - phi) - M * g * L1 * cos(theta_next);

            % NOTE: all other torques use current k values (Tac(k), T_balance(k), etc.)
            % muscle torque uses predicted MTF_next
            tau_other = F_k_up(k) + MTF_next * R + Tac(k) + T_balance(k) + tauGrav_next + tau_lim_next;

            Ieff = max(I, 1e-9);
            den  = 1 + (Ts / Ieff) * jointDamping(k);  % implicit viscous damping

            % implicit damping update for velocity
            dtheta_new = (dtheta(k) + (Ts / Ieff) * tau_other) / den;

            % semi-implicit position update
            theta_new  = theta(k) + Ts * dtheta_new;

            % --- MTU kinematics: update LMT using new velocity (consistent) ---
            dLMT_new = -R * dtheta_new;
            LMT_new  = LMT(k) + Ts * dLMT_new;

            % --- muscle-tendon force update (use predicted geometry, but new dLMT) ---
            % (this is still Euler-like, but inside iteration it becomes "implicit-ish")
            dFt_new  = k_t(k) * (dLMT_new + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k)));
            MTF_new  = MTF(k) + Ts * dFt_new;

            % --- relaxation to avoid oscillatory fixed-point (optional but helpful) ---
            relax = 0.6;  % 0.5~0.8 typical
            theta_next = relax * theta_next + (1 - relax) * theta_new;
            dtheta_next = relax * dtheta_next + (1 - relax) * dtheta_new;
            LMT_next   = relax * LMT_next   + (1 - relax) * LMT_new;
            MTF_next   = relax * MTF_next   + (1 - relax) * MTF_new;
        end

        % write back the converged-ish results
        theta(k + 1)  = theta_next;
        dtheta(k + 1) = dtheta_next;
        ddtheta(k)  = (dtheta(k + 1) - dtheta(k)) / Ts;

        dLMT(k + 1) = -R * dtheta(k + 1);
        LMT(k + 1)  = LMT_next;

        dFt(k)    = (MTF_next - MTF(k)) / Ts;
        MTF(k + 1)  = MTF_next;

        % store tau_lim and tauGrav based on final predicted state (optional)
        tau_lim(k) = local_tau_lim(theta(k), dtheta(k));  % or theta_next/dtheta_next if you prefer

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
    X(:, 7) = jointDamping;
    X(:, 8) = F_PE;
    X(:, 9) = Tac;
    X(:, 10) = T_balance;
    X(:, 11) = MTF(1:end - 1) - F_PE;

end

%% ===== Helper =====
% function [mLoad, I] = local_load_and_inertia(k, L1, L2, M)
%     if      k > 9060  && k < 12342, m = 5;
%     elseif  k > 22676 && k < 26140, m = 10;
%     elseif  k > 37973 && k < 40536, m = 15;
%     else,   m = 0;
%     end
%     mLoad = m;
%     I = m*L2^2 + (4*M*L1^2)/3;
% end

function [mLoad, I] = local_load_and_inertia(k, L1, L2, M)
    if      k > 6742  && k < 10348
        m = 5;
    elseif  k > 22885 && k < 25303
        m = 10;
    elseif  k > 37516 && k < 40124
        m = 15;
    elseif  k > 53975  && k < 57348
        m = 5;
    elseif  k > 67634 && k < 71055
        m = 10;
    elseif  k > 85027 && k < 87002
        m = 15;
    elseif  k > 99146  && k < 101117
        m = 5;
    elseif  k > 114510 && k < 117291
        m = 10;
    elseif  k > 128695 && k < 132744
        m = 15;
    else
        m = 0;
    end
    mLoad = m;
    I = m * L2^2 + (4 * M * L1^2) / 3;
end

function tau = local_tau_lim(theta, dtheta)
    theta_min = deg2rad(0);
    theta_max = deg2rad(95);

    a  = 30;
    k_wall = 1200;
    Fmax   = 900;
    c_wall = 0;

    H  = @(x) 0.5 * (1 + tanh(a * x));

    actUp   = H(theta - theta_max);
    actDown = H(theta_min - theta);

    penU = max(theta - theta_max, 0);
    penL = max(theta_min - theta, 0);

    tau_up   = -Fmax * tanh((k_wall / Fmax) * penU) * actUp;
    tau_down = +Fmax * tanh((k_wall / Fmax) * penL) * actDown;

    v_out_up   = max(dtheta,  0);
    v_out_down = max(-dtheta,  0);
    tau_damp   = -c_wall * (v_out_up * actUp - v_out_down * actDown);

    tau = tau_up + tau_down + tau_damp;
end
