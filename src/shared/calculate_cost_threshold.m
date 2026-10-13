function cost = calculate_cost_threshold(y_ref, y_model)
    % 1. 计算基础误差
    error = y_ref - y_model;

    % 2. 初始化权重向量，默认全为 1
    w = ones(size(y_ref));

    % 3. 【关键一步】找到所有 y_ref < 1 的点，把权重设大
    % 这里的 10 是惩罚倍数，你可以改成 20, 50 甚至 100
    w(y_ref < 1.4) = 10;

    % 4. 计算加权后的 Cost
    % 注意：点乘 (.*) 是为了让每个点对应自己的权重
    weighted_error = error .* w;

    cost = 10 * norm(weighted_error);
end
