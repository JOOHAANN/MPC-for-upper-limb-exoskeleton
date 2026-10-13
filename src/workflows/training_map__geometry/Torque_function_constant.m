% TORQUE_FUNCTION_CONSTANT Integrate the lumped muscle and exoskeleton dynamics.
% Use prepare_workflow to select matching code and input data.
function x = Torque_function_constant( ...
                                      Tac, theta_real, activation, Ts, M, L1, L2, phi, Torque0, LMT0)
    % TORQUE_FUNCTION_CONSTANT  (clean structured version of Program 1)
    % -------------------------------------------------------------------------
    % 输入:
    %   Tac          : 助力扭矩（列向量，N·m）
    %   theta_real   : 实际腰部角度 (rad)
    %   activation   : 肌肉激活 (0–1)   [在本版本中只用于长度]
    %   Ts           : 采样时间 (s)
    %   M            : 上半身质量 (kg)
    %   L1, L2       : 质心位置、外负载臂长 (m)
    %   phi          : 外负载相位角 (rad)
    %   Torque0      : 初始肌腱张力产生的扭矩（暂未使用，仅保留接口）
    %   LMT0         : 初始 MTU 长度（暂未使用，仅保留接口）
    %
    % 输出:
    %   x(:,1) = model_torque      (与原程序1相同)
    %   x(:,2) = theta_real
    %   x(:,3) = dtheta_real
    %   x(:,4) = ddtheta_real
    %   x(:,5) = F_k_up            (全 0，与原程序1一致)
    %   x(:,6) = T_front           (全 0)
    %   x(:,7) = activation        (原始激活)
    %   x(:,8) = F_PE              (全 0)
    %   x(:,9) = Tac
    %   x(:,10)= T_balance         (全 0)
    %   x(:,11)= MTF - F_PE        (此处 MTF = 0，保持原逻辑)
    % -------------------------------------------------------------------------

    %% ===== 常量 =====
    g  = 9.81;        % 重力加速度 [m/s^2]
    PercentageChange = 0.15;   %#ok<NASGU>  % 在本版本未使用, 仅保留
    c  = 60;          % 关节黏性阻尼(常数, 与原程序1一致)

    % 肌腱刚度多项式(原程序1中就有, 这里同样保留, 虽然本版本未使用)
    kt_spline = [ ...
                 6.07450471822638e+42  -1.17093401235721e+42  9.78861307802160e+40 ...
                 -4.68161016668128e+39   1.41111335336212e+38  -2.75023728531764e+36 ...
                 3.32535388381590e+34  -1.97635568280683e+32  -4.74852712911846e+29 ...
                 1.68812705982292e+28  -9.30681300424510e+25  -1.69227796153188e+23 ...
                 3.64048821058250e+21  -7.76789223923884e+18  -4.91004587706162e+16 ...
                 201470178531589        239829087287.116      -1720265398.11745 ...
                 -15739.3275509684       6708.90359529906      3.90612013684262]; %#ok<NASGU>

    %% ===== 预处理 & 预分配 =====
    activation = activation(:);     % 保证列向量
    Tac        = Tac(:);
    theta_real = theta_real(:);

    NumberOfData = length(activation);

    theta       = zeros(NumberOfData, 1);   %#ok<NASGU>
    dtheta      = zeros(NumberOfData, 1);   %#ok<NASGU>
    ddtheta     = zeros(NumberOfData, 1);   %#ok<NASGU>
    model_torque = zeros(NumberOfData, 1);

    LoM        = zeros(NumberOfData, 1);    %#ok<NASGU>
    alpha      = zeros(NumberOfData, 1);    %#ok<NASGU>
    dalpha     = zeros(NumberOfData, 1);    %#ok<NASGU>
    V_max      = zeros(NumberOfData, 1);    %#ok<NASGU>
    LMT        = zeros(NumberOfData, 1);    %#ok<NASGU>
    MA         = zeros(NumberOfData, 1);    %#ok<NASGU>
    LM         = zeros(NumberOfData, 1);    %#ok<NASGU>
    dLM        = zeros(NumberOfData, 1);    %#ok<NASGU>
    LT         = zeros(NumberOfData, 1);    %#ok<NASGU>
    dLMT       = zeros(NumberOfData, 1);    %#ok<NASGU>
    k_t        = zeros(NumberOfData, 1);    %#ok<NASGU>
    stiffness  = zeros(NumberOfData, 1);    %#ok<NASGU>
    F_L        = zeros(NumberOfData, 1);    %#ok<NASGU>
    F_PE       = zeros(NumberOfData, 1);
    t          = zeros(NumberOfData, 1);    %#ok<NASGU>
    F_k_up     = zeros(NumberOfData, 1);
    T_balance  = zeros(NumberOfData, 1);
    T_front    = zeros(NumberOfData, 1);
    MTF        = zeros(NumberOfData + 1, 1);  % 为兼容原输出
    dFt        = zeros(NumberOfData, 1);    %#ok<NASGU>
    strain     = zeros(NumberOfData, 1);    %#ok<NASGU>

    % 初始条件 (与原程序1保持同样赋值方式)
    alpha(1) = 0.23;          % pennation 初值
    LT(1)    = 0.05;          % 肌腱长度初值
    LMT(1)   = LMT0;
    dLMT(1)  = 0;
    dLM(1)   = 0;
    theta(1) = theta_real(1);

    % 角速度与角加速度(根据真实测量角度)
    dtheta_real  = [0; diff(theta_real) ./ Ts];
    ddtheta_real = [0; diff(dtheta_real) ./ Ts];

    %% ===== 主循环 =====
    for k = 1:NumberOfData - 1

        % --- 根据采样点选择负载与转动惯量 (与原程序1完全相同逻辑) ---
        [m, I] = local_load_and_inertia(k, L1, L2, M);

        % --- 关节净力矩 (完全照抄原 model_torque 公式) ---
        model_torque(k) = ddtheta_real(k) * I ...
                        + c * dtheta_real(k) ...
                        - Tac(k) ...
                        + m * g * L2 * cos(theta_real(k) - phi) ...
                        + M * g * L1 * cos(theta_real(k));

        % 原程序1在此处仅有被注释掉的 LMT 积分，故这里不做其它操作
    end

    % 最后一个点照旧保持为 0（与原程序1一致）
    model_torque(NumberOfData) = model_torque(NumberOfData - 1);

    %% ===== 组织输出，保持与程序1相同 =====
    x = zeros(NumberOfData, 11);
    x(:, 1)  = model_torque;
    x(:, 2)  = theta_real;
    x(:, 3)  = dtheta_real;
    x(:, 4)  = ddtheta_real;
    x(:, 5)  = F_k_up;
    x(:, 6)  = T_front;
    x(:, 7)  = activation;
    x(:, 8)  = F_PE;
    x(:, 9)  = Tac;
    x(:, 10) = T_balance;
    x(:, 11) = MTF(1:end - 1) - F_PE;   % 这里 MTF 为 0，与原程序1保持一致

end

%% ===== 局部函数：负载 & 转动惯量（完全照抄原逻辑） =====
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
