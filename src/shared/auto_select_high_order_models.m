function bestModel = auto_select_high_order_models(Baselines)
    % Baselines 是一个 200x4 的矩阵，其中 x 是 200x1，y 是 200x4

    x = Baselines(:, 1);  % 假设 Baselines 的第一列是 x
    y = Baselines(:, 2:5);  % 假设 Baselines 的第二到第五列是 y 的四个不同列

    % 定义可选模型和高阶多项式模型的最大阶数
    models = {'poly1', 'poly2', 'poly3', 'poly4', 'poly5', 'poly6', 'poly7', 'poly8', 'poly9', 'poly10'};
    numModels = length(models);
    numCols = size(y, 2);

    % 初始化存储最优模型的 cell 数组
    bestModel = cell(1, numCols);
    bestGof = zeros(1, numCols);

    % 对每一列的 y 进行拟合，并选择最优模型
    for i = 1:numCols
        bestGof(i) = -Inf;  % 初始化为负无穷大以便比较
        for j = 1:numModels
            % 进行拟合
            try
                [fitResult, gof] = fit(x, y(:, i), models{j});
                if gof.rsquare > bestGof(i)
                    bestGof(i) = gof.rsquare;
                    bestModel{i} = fitResult;
                end
            catch ME
                % 如果某个模型拟合失败，则继续尝试其他模型
                warning('Model %s failed for column %d: %s', models{j}, i, ME.message);
                continue
            end
        end

        % 显示最优模型的拟合参数
        fprintf('No %d column optimal model is: %s\n', i, formula(bestModel{i}));
    end
end
