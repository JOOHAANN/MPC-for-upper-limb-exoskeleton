function bestModel = auto_select_model(Baselines)
    % Baselines 是一个 200x4 的矩阵，其中 x 是 200x1，y 是 200x4

    x = Baselines(:, 1);  % 假设 Baselines 的第一列是 x
    y = Baselines(:, 2:5);  % 假设 Baselines 的第二到第五列是 y 的四个不同列

    % 平滑数据
    y = smoothdata(y, 'movmean', 5);

    % 定义可选模型
    models = {'poly1', 'poly2', 'exp1', 'exp2', 'gauss1', 'gauss2', 'power1', 'sin1', 'rat01', 'fourier1'};
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
                % 对某些模型设置初始参数
                if strcmp(models{j}, 'gauss1') || strcmp(models{j}, 'gauss2')
                    gaussModel = fittype('a*exp(-((x - b)^2) / (2*c^2)) + d', 'independent', 'x', 'coefficients', {'a', 'b', 'c', 'd'});
                    options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [1, mean(x), std(x), 0]);
                    [fitResult, gof] = fit(x, y(:, i), gaussModel, options);
                else
                    [fitResult, gof] = fit(x, y(:, i), models{j});
                end
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
        fprintf('第 %d 列的最优拟合模型为: %s\n', i, formula(bestModel{i}));
    end
end
