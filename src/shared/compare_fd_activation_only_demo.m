clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');
% FAST 10-second simulation: Activation-only Forward Dynamics (OpenSim-style)
% 唯一输入: activation(t)
% 输出: 关节角 theta(t) 与 肌肉力矩 tau(t)
% 模型与真实系统采用相同结构，参数略有偏差
% -> 仿真时长 10 秒

%% ---------- 唯一输入 ----------
u_fun = @(t) excitation_profile_long(t);  % 10s 激发模式

%% ---------- 参数 ----------
Ptrue = default_params_fast();
Pest  = default_params_fast();

% 参数轻微不匹配
Pest.Fmax   = 0.9  * Ptrue.Fmax;
Pest.r      = 1.05 * Ptrue.r;
Pest.Kpass  = 1.2  * Ptrue.Kpass;
Pest.B      = 0.8  * Ptrue.B;
Pest.Kgrav  = 0.95 * Ptrue.Kgrav;

% 激活动力学：稍大一点时间常数减少刚性
Ptrue.tau_act = 0.02;
Ptrue.tau_de = 0.05;
Pest.tau_act  = 0.022;
Pest.tau_de  = 0.055;

tEnd = 10.0;  % 仿真时长
x0   = [deg2rad(0); 0.0; 0.05];

% 加速用 ODE 设置
ode_opts = odeset('RelTol', 5e-4, 'AbsTol', 1e-6, 'MaxStep', 2e-2);

% 外扰（仅真实系统有）
tau_disturb_true = @(t) 0.5 * sin(2 * pi * 0.6 * t) .* (t > 3 & t < 8);

%% ---------- 仿真 ----------
[Tpl, Xpl, AuxPl] = simulate_fd_fast(Ptrue, u_fun, tEnd, x0, ode_opts, tau_disturb_true);
theta_true = Xpl(:, 1);
tau_true = AuxPl.tau_m;

rng(1);
theta_meas = theta_true + deg2rad(0.1) * randn(size(theta_true));
tau_meas   = tau_true   + 0.1 * randn(size(tau_true));

[Tm, Xm, AuxM] = simulate_fd_fast(Pest, u_fun, tEnd, x0, ode_opts, @(t)0);
theta_model  = Xm(:, 1);
tau_model    = AuxM.tau_m;

Tref = Tm;
theta_true_i = interp1(Tpl, theta_meas, Tref, 'linear', 'extrap');
tau_true_i   = interp1(Tpl, tau_meas, Tref, 'linear', 'extrap');

%% ---------- 误差 ----------
[rmse_theta, nrmse_theta] = err_metrics(theta_true_i, theta_model);
[rmse_tau,   nrmse_tau] = err_metrics(tau_true_i,   tau_model);

fprintf('\n=== 10 s simulation results ===\n');
fprintf('θ: RMSE = %.3f deg, NRMSE = %.2f %%\n', rad2deg(rmse_theta), nrmse_theta * 100);
fprintf('τ: RMSE = %.3f N·m, NRMSE = %.2f %%\n', rmse_tau, nrmse_tau * 100);

%% ---------- 绘图 ----------
figure('Color', 'w', 'Position', [80 80 1100 720]);

subplot(3, 1, 1);
plot(Tref, rad2deg(theta_true_i), 'k', 'LineWidth', 1.4);
hold on;
plot(Tref, rad2deg(theta_model), '--', 'LineWidth', 1.4);
grid on;
ylabel('\theta (deg)');
legend('True', 'Model');
title('Joint angle');

subplot(3, 1, 2);
plot(Tref, tau_true_i, 'k', 'LineWidth', 1.4);
hold on;
plot(Tref, tau_model, '--', 'LineWidth', 1.4);
grid on;
ylabel('\tau (N·m)');
legend('True', 'Model');
title('Torque');

subplot(3, 1, 3);
plot(Tref, rad2deg(theta_model - theta_true_i), 'LineWidth', 1.0);
hold on;
plot(Tref, tau_model - tau_true_i, '--', 'LineWidth', 1.0);
grid on;
xlabel('t (s)');
ylabel('Error');
legend('\Delta\theta (deg)', '\Delta\tau (N·m)');
title('Model − True');

sgtitle('10 s Activation-only Forward Dynamics (Fast OpenSim-style)');

%% ======= 支撑函数 =======
function [T, X, Aux] = simulate_fd_fast(P, u_fun, tEnd, x0, ode_opts, tau_disturb_fun)
    if nargin < 6
        tau_disturb_fun = @(t)0;
    end
    rhs = @(t, x) fd_rhs_fast(t, x, P, u_fun, tau_disturb_fun);
    [T, X] = ode23tb(rhs, [0 tEnd], x0, ode_opts);
    n  = numel(T);
    tau_m = zeros(n, 1);
    tau_g = tau_m;
    tau_p = tau_m;
    tau_d = tau_m;
    for k = 1:n
        [~, aux] = rhs(T(k), X(k, :).');
        tau_m(k) = aux.tau_m;
        tau_g(k) = aux.tau_g;
        tau_p(k) = aux.tau_p;
        tau_d(k) = aux.tau_d;
    end
    Aux = struct('tau_m', tau_m, 'tau_g', tau_g, 'tau_p', tau_p, 'tau_d', tau_d);
end

function [dx, aux] = fd_rhs_fast(t, x, P, u_fun, tau_disturb_fun)
    q  = x(1);
    dq = x(2);
    a = x(3);

    % 唯一输入：activation/excitation
    u = u_fun(t);
    H = double(u > a);
    tau_a = P.tau_de + (P.tau_act - P.tau_de) * H;
    a = min(max(a, 0.0), 1.0);                 % ← 夹紧激活
    da = (u - a) / tau_a;

    % 肌肉几何
    Lmt = P.Lmt0 - P.r * (q - P.q0);
    Vmt =           -P.r * dq;

    % 刚性腱近似 + 纤维长度夹紧（防爆）
    Lf  = Lmt - P.Lslack;
    Lf  = min(max(Lf, 0.20 * P.Lopt), 2.5 * P.Lopt);   % ← 新增上界
    Vf  = Vmt;

    % 力-长/力-速
    Fl = exp(-((Lf - P.Lopt)^2) / (2 * (P.FLw^2)));
    v  = Vf / P.vmax;
    A = 0.25;
    Fv = (1 - v) / (1 + v / A);
    Fv = min(max(Fv, 0.05), 1.8);

    % 被动纤维（稳定化：expm1 + 上界）
    xpe = Lf / P.Lopt - 1.0;
    if xpe > P.Fpe_s
        arg = P.Fpe_k * (xpe - P.Fpe_s);
        arg = min(arg, 8.0);                      % ← 饱和，防溢出
        num = expm1(arg);
        den = expm1(P.Fpe_k * (1 - P.Fpe_s));
        Fpe = (num / den) * 0.4;
    else
        Fpe = 0.0;
    end

    Fm  = P.Fmax * (a * Fl * Fv + Fpe);
    tau_m = P.r * Fm;

    % 重力/被动/软限位
    tau_g = P.Kgrav * sin(q);
    tau_p = P.Kpass * (q - P.q0) + P.B * dq;

    % 软限位（可选：若你已经很稳，可去掉）
    q_min = deg2rad(-60);
    q_max = deg2rad(60);
    k_lim = 30;
    b_lim = 2.0;
    tau_lim = k_lim * (tanh(10 * (q - q_max)) - tanh(10 * (q - q_min))) ...
              + b_lim * max(min(dq, 5), -5);

    % 外扰
    tau_d = tau_disturb_fun(t);

    % 动力学
    qdd = (tau_m - tau_g - tau_p - tau_lim + tau_d) / P.I;

    dx = [dq; qdd; da];
    if nargout > 1
        aux = struct('tau_m', tau_m, 'tau_g', tau_g, 'tau_p', tau_p, ...
                     'tau_d', tau_d, 'tau_lim', tau_lim);
    end
end

function P = default_params_fast()
    P.I = 1.0;
    P.q0 = 0;
    P.B = 1.0;
    P.Kpass = 5;
    P.Kgrav = 50;
    P.r = 0.04;
    P.Lopt = 0.1;
    P.Lslack = 0.2;
    P.Lmt0 = P.Lopt + P.Lslack;
    P.Fmax = 1500;
    P.vmax = 10 * P.Lopt;
    P.FLw = 0.45 * P.Lopt;
    P.Fpe_k = 3;
    P.Fpe_s = 0.6;
    P.tau_act = 0.02;
    P.tau_de = 0.05;
end

function u = excitation_profile_long(t)
    % 10 秒内的平滑激发
    if t < 2
        u = 0.1 + 0.05 * sin(pi * t / 2);
    elseif t < 5
        u = 0.4 + 0.05 * sin(2 * pi * (t - 2) / 3);
    elseif t < 8
        u = 0.25 + 0.03 * sin(2 * pi * (t - 5) / 3);
    else
        u = 0.15 + 0.02 * sin(2 * pi * (t - 8) / 2);
    end
end

function [rmse, nrmse] = err_metrics(y_true, y_hat)
    rmse = sqrt(mean((y_hat - y_true).^2));
    range = max(y_true) - min(y_true);
    nrmse = rmse / max(range, 1e-12);
end
