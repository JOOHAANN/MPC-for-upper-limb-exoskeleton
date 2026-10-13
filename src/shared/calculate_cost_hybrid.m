function total_cost = calculate_cost_hybrid(y_ref, y_model)

    % 强制成列向量，避免 corrcoef/广播的形状坑
    y_ref   = y_ref(:);
    y_model = y_model(:);

    % 长度检查（不一致就直接惩罚）
    if numel(y_ref) ~= numel(y_model) || isempty(y_ref)
        total_cost = 1e12;
        return
    end

    % 实数检查（出现复数就惩罚）
    if ~isreal(y_ref) || ~isreal(y_model)
        total_cost = 1e12;
        return
    end

    % 有限性检查（NaN/Inf 就惩罚）
    if any(~isfinite(y_ref)) || any(~isfinite(y_model))
        total_cost = 1e12;
        return
    end

    % ===== 1) Correlation cost =====
    R_matrix = corrcoef(y_ref, y_model);
    R = R_matrix(1, 2);
    if ~isfinite(R)   % 包括 NaN/Inf
        R = 0;
    end
    cost_correlation = 1 - R;

    % ===== 2) Weighted error cost =====
    error = y_ref - y_model;
    w = ones(size(y_ref));
    w(y_ref < 1) = 10;
    cost_error = norm(error .* w);

    % ===== 3) Total =====
    lambda = 50000;
    total_cost = 10 * cost_error + lambda * cost_correlation;

    % 最终兜底：确保返回“实数标量”
    if ~isscalar(total_cost) || ~isreal(total_cost) || ~isfinite(total_cost)
        total_cost = 1e12;
    end
end
