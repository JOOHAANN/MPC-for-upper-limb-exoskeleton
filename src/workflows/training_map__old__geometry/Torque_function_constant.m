function x = Torque_function_constant(Tac, theta_real, activation, Ts, M, L1, L2, phi, Torque0, LMT0)
    % TORQUE_FUNCTION_CONSTANT Integrate the lumped muscle and exoskeleton dynamics.
    % Use prepare_workflow to select matching code and input data.
    % TORQUE_FUNCTION_CONSTANT
    % 计算基于测得角度与外部负载的模型力矩及相关通道
    %
    % 输入：
    %   Tac          : 外骨骼/电机输出力矩（Nx1 / 1xN）
    %   theta_real   : 实测关节角（rad）（Nx1 / 1xN）
    %   activation   : 肌激活（Nx1 / 1xN）— 这里仅作为通道输出，不参与计算
    %   Ts           : 采样周期（s）
    %   M            : 上躯干等等效质量（kg）
    %   L1, L2       : 几何长度参数（m）
    %   phi          : 前方载荷与躯干的角度偏置（rad）
    %   Torque0      : 预留参数（未使用）
    %   LMT0         : 初始肌腱-肌腹长度（m）（保留）
    %
    % 输出：
    %   x : [N x 11] 矩阵，各列含义见文末

    %% --------- 统一输入为列向量 & 基本检查 ---------
    Tac        = Tac(:);
    theta_real = theta_real(:);
    activation = activation(:);

    N = numel(theta_real);

    assert(numel(Tac) == N,        'Tac 长度需与 theta_real 一致。');
    assert(numel(activation) == N, 'activation 长度需与 theta_real 一致。');
    assert(Ts > 0,                 'Ts 必须为正。');

    g = 9.81;

    %% --------- 预分配（全部列向量）---------
    theta        = zeros(N, 1); %#ok<NASGU>  % 如需后续积分可启用
    dtheta       = zeros(N, 1); %#ok<NASGU>
    ddtheta      = zeros(N, 1); %#ok<NASGU>
    dtheta_real  = zeros(N, 1);
    ddtheta_real = zeros(N, 1);
    model_torque = zeros(N, 1);

    % 下方通道目前未实际使用，保持与原版兼容
    LoM   = zeros(N, 1); %#ok<NASGU>
    alpha = zeros(N, 1);
    dalpha = zeros(N, 1); %#ok<NASGU>
    V_max = zeros(N, 1); %#ok<NASGU>
    LMT   = zeros(N, 1);
    MA    = zeros(N, 1); %#ok<NASGU>
    LM    = zeros(N, 1); %#ok<NASGU>
    dLM   = zeros(N, 1); %#ok<NASGU>
    LT    = zeros(N, 1);
    dLMT  = zeros(N, 1); %#ok<NASGU>
    k_t   = zeros(N, 1); %#ok<NASGU>
    stiffness = zeros(N, 1); %#ok<NASGU>
    F_L   = zeros(N, 1); %#ok<NASGU>
    F_PE  = zeros(N, 1);
    t     = (0:N - 1).' * Ts; %#ok<NASGU> % 时间轴
    F_k_up   = zeros(N, 1);
    F_k_down = zeros(N, 1);
    T_balance = zeros(N, 1);
    T_front  = zeros(N, 1);

    MTF  = zeros(N + 1, 1);      % 与原版保持 N+1
    dFt  = zeros(N, 1); %#ok<NASGU>
    strain = zeros(N, 1); %#ok<NASGU>

    % 预设值（保持与你原始代码兼容）
    alpha(1) = 0.23;
    LT(1)    = 0.05;
    LMT(1)   = LMT0;

    %% --------- 数值微分（鲁棒边界处理）---------
    % 采用差分并与梯度做一次融合，避免首末端过尖
    dtheta_tmp       = [0; diff(theta_real)] / Ts;
    ddtheta_tmp      = [0; diff(dtheta_tmp)] / Ts;
    dtheta_real      = dtheta_tmp;
    ddtheta_real     = ddtheta_tmp;

    %% --------- 阻尼与模型常数（可按需要暴露为参数）---------
    c_damp = 10;          % 阻尼
    % a = 0.92;           % 原代码未用，移除以减少警告

    %% --------- 主循环：计算 model_torque ---------
    for k = 1:N
        % 当前时刻的外部负载与惯量（已修复为“按时间”而非固定下标）
        [m_ext, I] = local_load_and_inertia(k, L1, L2, M, N);

        % 力矩平衡（你的原方程保持不变）
        %  I*ddtheta + c*dtheta - Tac + m*g*L2*cos(theta-phi) + M*g*L1*cos(theta) = model_torque
        model_torque(k) = ddtheta_real(k) * I + c_damp * dtheta_real(k) ...
                        - Tac(k) ...
                        + m_ext * g * L2 * cos(theta_real(k) - phi) ...
                        + M * g * L1 * cos(theta_real(k));
    end

    % 平滑（窗口自动限制不超过 N，避免窗口大于序列长度报错）
    if N >= 5
        win = min(300, max(5, 2 * ceil(0.05 * N / 2) + 1));  % 经验：5%长度的奇数窗
        model_torque = smoothdata(model_torque, "gaussian", win);
    end

    % 兼容输出：LM_LoM 未在后续使用，仅保留写法避免下游依赖断裂
    LM_LoM = zeros(N, 1); %#ok<NASGU>

    %% --------- 汇总输出 ---------
    x = zeros(N, 11);
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
    x(:, 11) = MTF(1:end - 1) - F_PE;  % 与原逻辑一致

end % ===== 主函数结束 =====

%% ================== 本地函数 ==================
% function [m, I] = local_load_and_inertia(k, L1, L2, M, N)
%
%     persistent mSmooth
%
%     if isempty(mSmooth)
%         % ===== 第一次调用时预计算整段质量曲线 =====
%         Nstep = N;              % 最大 k，按实际改
%
%         mStep = zeros(Nstep,1);
%         mStep(9048:11971)  = 5;
%         mStep(22700:26083) = 10;
%         mStep(38023:40518) = 15;
%
%         Fs   = 100;                 % 采样频率
%         Fc   = 2;                   % 截止频率
%         [b,a] = butter(2, Fc/(Fs/2), 'low');
%         mSmooth = filtfilt(b,a,mStep);
%         mSmooth(mSmooth<0) = 0;
%     end
%
%     % ===== 每次调用时只查当前 k =====
%     k = max(1, min(k, length(mSmooth)));
%     m = mSmooth(k);
%
%     I = m*L2^2 + (4*M*L1^2)/3;
% end
function [m, I] = local_load_and_inertia(k, L1, L2, M, N)
    % 在给定的矩形区间上加“倒圆角”——圆角过渡允许越过蓝线（区间）外侧
    % 红线不会被夹在区间内，而是在开始点前先起、结束点后再落。

    persistent mSmooth_cached N_cached

    if isempty(mSmooth_cached) || isempty(N_cached) || N_cached ~= N
        % ===== 你给的“蓝线”矩形区间与幅值 A =====
        pulses = [
                  9048   11971    5
                  22700  26083   10
                  38023  40518   15
                 ];

        % ===== 圆角 & 多项式参数 =====
        edgeFrac  = 0.10;  % 圆角半径占区间长度的比例（两边各 R），大=更圆润、更越界
        polyOrder = 7;     % 3/5/7（阶越高越丝滑；7 为 C^3）

        mSmooth_cached = build_rounded_box(N, pulses, edgeFrac, polyOrder);
        N_cached = N;
    end

    % ===== 查询当前索引 =====
    k = max(1, min(k, N_cached));
    m = mSmooth_cached(k);

    % ===== 转动惯量 =====
    I = m * L2^2 + (4 * M * L1^2) / 3;
end

% ----------------- 辅助实现 -----------------
function y = build_rounded_box(N, pulses, edgeFrac, polyOrder)
    y = zeros(N, 1);

    for i = 1:size(pulses, 1)
        s = max(1, min(N, pulses(i, 1)));   % 区间开始（蓝线）
        e = max(1, min(N, pulses(i, 2)));   % 区间结束（蓝线）
        if e < s
            [s, e] = deal(e, s);
        end
        A = pulses(i, 3);

        L = e - s + 1;                     % 区间宽度
        if L <= 2
            continue
        end

        % 圆角半径（越界宽度）：两侧各 R
        R = max(1, round(edgeFrac * L));

        % —— 左圆角：从 s-R 到 s+R，用多项式把 0 -> A
        L1 = max(1, s - R);
        L2 = min(N, s + R);
        if L2 >= L1
            n  = L2 - L1 + 1;
            t  = linspace01(n);            % [0,1]
            sCurve = poly_smoothstep(t, polyOrder);
            seg = A * sCurve;
            y(L1:L2) = max(y(L1:L2), seg);
        end

        % —— 右圆角：从 e-R 到 e+R，用多项式把 A -> 0
        R1 = max(1, e - R);
        R2 = min(N, e + R);
        if R2 >= R1
            n  = R2 - R1 + 1;
            t  = linspace01(n);
            sCurve = 1 - poly_smoothstep(t, polyOrder);
            seg = A * sCurve;
            y(R1:R2) = max(y(R1:R2), seg);
        end

        % —— 中间平台：在 [s+R, e-R] 之间保持 A
        coreStart = s + R;
        coreEnd   = e - R;
        if coreStart <= coreEnd
            coreStart = max(1, min(N, coreStart));
            coreEnd   = max(1, min(N, coreEnd));
            y(coreStart:coreEnd) = max(y(coreStart:coreEnd), A);
        end
    end
end

function s = poly_smoothstep(t, order)
    % 多项式过渡（端点及若干阶导为 0）
    % order=3:  3t^2 - 2t^3                    (C^1)
    % order=5:  6t^5 - 15t^4 + 10t^3           (C^2)
    % order=7:  35t^4 - 84t^5 + 70t^6 - 20t^7  (C^3)
    switch order
        case 3
            s = t.^2 .* (3 - 2 * t);
        case 5
            s = t.^3 .* (10 - 15 * t + 6 * t.^2);
        case 7
            s = t.^4 .* (35 - 84 * t + 70 * t.^2 - 20 * t.^3);
        otherwise
            s = t.^4 .* (35 - 84 * t + 70 * t.^2 - 20 * t.^3);
    end
end

function t = linspace01(n)
    % 生成 [0,1] 间等距列；n=1 时返回 1（避免除零）
    if n <= 1
        t = 1;
    else
        t = (0:n - 1).' / (n - 1);
    end
end
