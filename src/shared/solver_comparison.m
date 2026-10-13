function solver_comparison
    clc;
    clear;
    close all;

    %% 1. 系统配置
    % 模拟一个刚性很大的肌肉-肌腱系统
    % 方程: m*x'' + c*x' + k*x^3 = F_drive * sin(...)
    p.m = 1.0;
    p.c = 5.0;       % 阻尼 (适中)
    p.k = 8000;      % 刚度很大 (这是关键！刚度越大，nIt=1死得越快)

    % 激励
    p.Amp = 500;
    p.Freq = 1.0;    % 1Hz 正弦波

    T_end = 5.0;     % 仿真5秒
    Ts = 0.01;       % 10ms 步长 (对于刚性系统，这个步长其实很危险)
    t_fixed = 0:Ts:T_end;

    y0 = [0; 0];     % 初始静止

    %% 2. 运行求解器

    % --- A. 标准答案 (ODE15s) ---
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t_15s, y_15s] = ode15s(@(t, y) sys_dyn(t, y, p), [0 T_end], y0, opts);

    % --- B. 你的方案：只迭代 1 次 (nIt = 1) ---
    % 这种情况下，力 F_spring 只能基于“当前时刻”计算 (Explicit Stiffness)
    y_iter1 = zeros(2, length(t_fixed));
    y_iter1(:, 1) = y0;

    tic;
    for k = 1:length(t_fixed) - 1
        x_curr = y_iter1(1, k);
        v_curr = y_iter1(2, k);
        t_now  = t_fixed(k);

        F_ext = p.Amp * sin(2 * pi * p.Freq * t_now);

        % [关键差异点] nIt=1 时，必须用 x_curr 算力
        F_spring = p.k * (x_curr)^3;

        % 隐式阻尼更新 (依然保留)
        v_next = (v_curr + (Ts / p.m) * (F_ext - F_spring)) / (1 + Ts * p.c / p.m);
        x_next = x_curr + Ts * v_next;

        y_iter1(1, k + 1) = x_next;
        y_iter1(2, k + 1) = v_next;
    end
    time_1 = toc;

    % --- C. 你的方案：迭代 3 次 (nIt = 3) ---
    % 这种情况下，力 F_spring 是基于“预测时刻”计算的 (Implicit Stiffness)
    y_iter3 = zeros(2, length(t_fixed));
    y_iter3(:, 1) = y0;

    tic;
    for k = 1:length(t_fixed) - 1
        x_curr = y_iter3(1, k);
        v_curr = y_iter3(2, k);
        t_now  = t_fixed(k);

        F_ext = p.Amp * sin(2 * pi * p.Freq * t_now);

        % 预测初始值
        x_pred = x_curr;
        v_pred = v_curr;

        % 内循环
        for it = 1:3
            % [关键差异点] 使用预测位置 x_pred 算力
            F_spring = p.k * (x_pred)^3;

            v_new = (v_curr + (Ts / p.m) * (F_ext - F_spring)) / (1 + Ts * p.c / p.m);
            x_new = x_curr + Ts * v_new;

            % 松弛
            relax = 0.6;
            x_pred = relax * x_pred + (1 - relax) * x_new;
            v_pred = relax * v_pred + (1 - relax) * v_new;
        end

        y_iter3(1, k + 1) = x_pred;
        y_iter3(2, k + 1) = v_pred;
    end
    time_3 = toc;

    %% 3. 结果绘图
    figure('Color', 'w', 'Position', [100, 100, 1200, 600]);

    % 子图1：位置对比
    subplot(2, 1, 1);
    hold on;
    grid on;
    plot(t_15s, y_15s(:, 1), 'k', 'LineWidth', 2, 'DisplayName', 'ODE15s (真值)');
    plot(t_fixed, y_iter3(1, :), 'b-o', 'MarkerSize', 4, 'DisplayName', 'nIt=3 (你的原方案)');
    plot(t_fixed, y_iter1(1, :), 'r.-', 'LineWidth', 1, 'DisplayName', 'nIt=1 (只迭代一次)');

    ylabel('位置 x (m)');
    title(['刚度 k=' num2str(p.k) ', 步长 Ts=' num2str(Ts) 's']);
    legend('Location', 'best');

    % 子图2：误差对比
    subplot(2, 1, 2);
    hold on;
    grid on;
    % 对齐数据
    y_true_interp = interp1(t_15s, y_15s(:, 1), t_fixed);

    err_1 = abs(y_iter1(1, :) - y_true_interp);
    err_3 = abs(y_iter3(1, :) - y_true_interp);

    plot(t_fixed, err_1, 'r-', 'LineWidth', 1.5, 'DisplayName', 'nIt=1 的误差');
    plot(t_fixed, err_3, 'b-', 'LineWidth', 1.5, 'DisplayName', 'nIt=3 的误差');

    ylabel('绝对误差 (m)');
    xlabel('时间 (s)');
    title('误差随时间的变化');
    legend;

    fprintf('nIt=1 耗时: %.5f s\n', time_1);
    fprintf('nIt=3 耗时: %.5f s\n', time_3);
end

function dydt = sys_dyn(t, y, p)
    x = y(1);
    F_ext = p.Amp * sin(2 * pi * p.Freq * t);
    dydt = [y(2); (F_ext - p.c * y(2) - p.k * x^3) / p.m];
end
