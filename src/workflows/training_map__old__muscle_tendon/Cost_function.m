function Cost = Cost_function(theta_real_1, Torque_real_1, Smooth_activation_1, Smooth_activation_front, damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, M, L1, L2, phi, theta0, Torque0, LMT0)
    % COST_FUNCTION Evaluate the calibration objective for this model variant.
    % Use prepare_workflow to select matching code and input data.

    % 一阶、二阶导（目前没在 cost 里用到，但保留以防 Torque_function_constant 里面会用）
    dtheta_real_1  = [0; diff(theta_real_1) ./ Ts];
    ddtheta_real_1 = [0; diff(dtheta_real_1) ./ Ts];

    % 这里先不给主动肌力（全 0），保持你原来的设定
    Tac = zeros(length(Smooth_activation_1), 1)';

    % 仿真得到输出：X1(:,1) = torque_hat, X1(:,2) = theta_hat
    X1 = Torque_function_constant( ...
                                  Tac, theta_real_1, Torque_real_1, ...
                                  Smooth_activation_1, Smooth_activation_front, ...
                                  damping, alpha0, Ts, R, LoM_nominal, LsT, FoM, ...
                                  M, L1, L2, phi, theta0, Torque0, LMT0);

    k = 1:length(Torque_real_1);
    valid_idx = (k > 7757 & k < 13853) | ...
                (k > 20917 & k < 27659) | ...
                (k > 34606 & k < 38513);
    % 另一种划分方式留着：
    % valid_idx = (k > 7757 & k < 10319) | (k > 11144 & k < 13828) | ...
    %             (k > 20956 & k < 23873) | (k > 24955 & k < 26758) | ...
    %             (k > 36908 & k < 39141) | (k > 39753 & k < 41240);

    %% ---------- 1. Torque：只看误差（RMSE） ----------
    yT_true = Torque_real_1(valid_idx);
    yT_hat  = X1(valid_idx, 1);

    eT   = yT_hat - yT_true;
    rmse_T = sqrt(mean(eT.^2));      % torque 的均方根误差

    % 可以保留你原来的 0.02 权重，也可以再微调
    w_torque = 0.02;
    Cost_torque = w_torque * rmse_T;

    %% ---------- 2. Theta：误差 + 形状（相关性） ----------
    yTh_true = theta_real_1(valid_idx);
    yTh_hat  = X1(valid_idx, 2);

    eTh = yTh_hat - yTh_true;

    % ① 误差项：归一化 RMSE（避免幅值太大自动 cost 大）
    rmse_Th = sqrt(mean(eTh.^2));
    rmse_Th_norm = rmse_Th / (std(yTh_true) + eps);   % 用 std 归一化

    % ② 形状项：相关系数（越接近 1 越好）
    C = corrcoef(yTh_hat, yTh_true);      % 2x2 矩阵
    if numel(C) >= 4
        r = C(1, 2);                       % 相关系数 r ∈ [-1, 1]
    else
        r = 0;                            % 万一数据点极少时，避免出错
    end

    % ③ 权重组合（你可以根据感觉再调）
    w_rmse = 0.7;      % 更看重数值误差
    w_corr = 0.3;      % 也看重形状相似度

    Cost_theta = w_rmse * rmse_Th_norm + w_corr * (1 - r);
    % 如果想让 Cost_theta 在数值上更大一点，方便和 torque 权重平衡，也可以再乘一个系数：
    % scale_theta = 100;
    % Cost_theta = scale_theta * Cost_theta;

    %% ---------- 总 cost ----------
    Cost = Cost_theta + Cost_torque;

    %% ---------- 数值保护：防止 NaN / Inf ----------
    if isnan(Cost) || isinf(Cost)
        Cost = 1e9;
    end

end
