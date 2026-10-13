function demo_ekf_rts_torque_measurement()
    % EKF + RTS smoother demo (torque as measurement, theta as state)
    % model: tau = I*ddtheta + c*dtheta - Tac + m*g*L2*cos(theta-phi) + M*g*L1*cos(theta)

    %% 0) Settings & constants
    rng(0);
    N  = 400;                 % steps
    dt = 0.01;                % sample time
    t  = (0:N - 1)' * dt;

    c   = 10;
    m   = 5.0;
    M   = 22.9819;
    g   = 9.81;
    L1  = 0.345;
    L2  = 0.756835;
    phi = 0.950692;
    I   = m * L2^2 + (4 * M * L1^2) / 3;
    Tac = 0;                  % no actuation

    % angle limits (user gave 1.6, 0.5 in reverse order)
    tmp = sort([1.6, 0.5]);   % ensure lb <= ub
    theta_lb = tmp(1);
    theta_ub = tmp(2);

    %% 1) Make W-shaped truth (only for validation)
    [theta_true, dtheta_true, ddtheta_true] = make_W_theta(t, 0.18, theta_lb, theta_ub, 1.0);
    tau_true = I * ddtheta_true + c * dtheta_true ...
             + m * g * L2 * cos(theta_true - phi) + M * g * L1 * cos(theta_true) - Tac;
    tau = tau_true + 2.0 * randn(N, 1);    % noisy torque measurement

    %% 2) EKF forward pass
    % state x = [theta; dtheta]
    x_est = [(theta_lb + theta_ub) / 2; 0];
    P     = diag([0.05^2, 1.0^2]);

    Q = diag([1e-5, 1e-3]);   % process noise
    R = 4.0;                  % torque measurement noise variance

    x_fwd  = zeros(2, N);
    P_fwd  = zeros(2, 2, N);
    x_pred = zeros(2, N);
    P_pred = zeros(2, 2, N);

    for k = 1:N
        if k == 1
            xk = x_est;
            Pk = P;
            tau_prev = tau(1);
        else
            xk = x_fwd(:, k - 1);
            Pk = P_fwd(:, :, k - 1);
            tau_prev = tau(k - 1);
        end

        % ---- Predict (Euler integration of dynamics) ----
        th  = xk(1);
        dth = xk(2);
        ddth = (tau_prev + Tac - c * dth ...
                - m * g * L2 * cos(th - phi) - M * g * L1 * cos(th)) / I;

        xk_pred = [th + dth * dt
                   dth + ddth * dt];

        dddth_dth = (m * g * L2 * sin(th - phi) + M * g * L1 * sin(th)) / I;
        dddth_dd  = (-c) / I;
        F = [1,           dt
             dddth_dth * dt, 1 + dddth_dd * dt];

        Pk_pred = F * Pk * F' + Q;

        % ---- Update with torque measurement ----
        yk = tau(k);
        thp  = xk_pred(1);
        dthp = xk_pred(2);
        h    = c * dthp + m * g * L2 * cos(thp - phi) + M * g * L1 * cos(thp) - Tac;
        dh_dth = -m * g * L2 * sin(thp - phi) - M * g * L1 * sin(thp);
        dh_dd  =  c;
        H = [dh_dth, dh_dd];

        S = H * Pk_pred * H' + R;
        K = (Pk_pred * H') / S;
        xk_upd = xk_pred + K * (yk - h);
        Pk_upd = (eye(2) - K * H) * Pk_pred;

        % angle bounds projection
        xk_upd(1) = min(max(xk_upd(1), theta_lb), theta_ub);

        % store
        x_pred(:, k) = xk_pred;
        P_pred(:, :, k) = Pk_pred;
        x_fwd(:, k)  = xk_upd;
        P_fwd(:, :, k)  = Pk_upd;
    end

    %% 3) RTS smoother (backward)
    x_smooth = x_fwd;
    P_smooth = P_fwd;
    for k = N - 1:-1:1
        th  = x_fwd(1, k);
        dth = x_fwd(2, k);
        dddth_dth = (m * g * L2 * sin(th - phi) + M * g * L1 * sin(th)) / I;
        dddth_dd  = (-c) / I;
        Fk = [1,           dt
              dddth_dth * dt, 1 + dddth_dd * dt];
        Ck = P_fwd(:, :, k) * Fk' / P_pred(:, :, k + 1);
        x_smooth(:, k)   = x_fwd(:, k) + Ck * (x_smooth(:, k + 1) - x_pred(:, k + 1));
        P_smooth(:, :, k) = P_fwd(:, :, k) + Ck * (P_smooth(:, :, k + 1) - P_pred(:, :, k + 1)) * Ck.';
        x_smooth(1, k) = min(max(x_smooth(1, k), theta_lb), theta_ub);
    end

    %% 4) Report & plots
    theta_ekf = x_fwd(1, :)';
    dtheta_ekf = x_fwd(2, :)';
    theta_rts = x_smooth(1, :)';
    dtheta_rts = x_smooth(2, :)';

    rmse_ekf = sqrt(mean((theta_ekf - theta_true).^2));
    rmse_rts = sqrt(mean((theta_rts - theta_true).^2));
    fprintf('RMSE(theta): EKF=%.4f  RTS=%.4f (rad)\n', rmse_ekf, rmse_rts);

    figure('Color', 'w', 'Position', [100 80 1100 820]);
    subplot(3, 1, 1);
    plot(t, theta_true, 'k-', 'LineWidth', 1.6);
    hold on;
    plot(t, theta_ekf, 'r--', 'LineWidth', 1.3);
    plot(t, theta_rts, 'b-.', 'LineWidth', 1.3);
    yline(theta_lb, 'b:', 'lb');
    yline(theta_ub, 'b:', 'ub');
    xlabel('Time (s)');
    ylabel('\theta (rad)');
    legend('\theta true', 'EKF', 'RTS smooth', 'Location', 'best');
    title('\theta estimation (torque-only, Tac=0)');

    subplot(3, 1, 2);
    plot(t, tau, 'k-', 'LineWidth', 1.0);
    hold on;
    plot(t, c * dtheta_rts + m * g * L2 * cos(theta_rts - phi) + M * g * L1 * cos(theta_rts) - Tac, ...
         'r--', 'LineWidth', 1.0);
    xlabel('Time (s)');
    ylabel('Torque (N·m)');
    legend('measured \tau', 'h(x_{RTS})', 'Location', 'best');
    title('Torque measurement vs h(x)');

    subplot(3, 1, 3);
    plot(t, dtheta_true, 'k--', 'LineWidth', 1.0);
    hold on;
    plot(t, dtheta_ekf, 'r-', 'LineWidth', 1.0);
    plot(t, dtheta_rts, 'b-', 'LineWidth', 1.0);
    xlabel('Time (s)');
    ylabel('d\theta (rad/s)');
    legend('d\theta true', 'EKF', 'RTS', 'Location', 'best');
    title('Angular velocity');

    sgtitle('EKF + RTS with torque measurement');

end

% ===== helpers =====
function [theta_true, dtheta_true, ddtheta_true] = make_W_theta(t, w, theta_lb, theta_ub, ascl)
    N  = numel(t);
    dt = t(2) - t(1);
    t1 = t(round(N * 0.30));
    t2 = t(round(N * 0.70));
    g1 = exp(-0.5 * ((t - t1) / w).^2);
    g2 = exp(-0.5 * ((t - t2) / w).^2);
    r  = g1 + g2;
    r = (r - min(r)) / (max(r) - min(r) + eps);
    theta_peak   = theta_ub - 0.05;
    theta_valley = theta_lb + 0.05;
    amp = (theta_peak - theta_valley) * ascl;
    theta_true = theta_peak - amp * r;
    theta_true = min(max(theta_true, theta_lb), theta_ub);
    dtheta_true  = dtheta_from_theta(theta_true, dt);
    ddtheta_true = D2_operator(N, dt) * theta_true;
end

function dthN = dtheta_from_theta(th, dt)
    N = numel(th);
    dthN = zeros(N, 1);
    dthN(1)     = (th(2)   - th(1))   / dt;
    dthN(2:N - 1) = (th(3:N) - th(1:N - 2)) / (2 * dt);
    dthN(N)     = (th(N)   - th(N - 1)) / dt;
end

function D2 = D2_operator(N, dt)
    e  = ones(N, 1);
    D2 = spdiags([e -2 * e e], -1:1, N, N) / (dt^2);
end
