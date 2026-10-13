function bestModel = excitation_primitives_double_gauss(Baselines)

    %% %%%%%%%%%%%%%%%%%
    x = Baselines(:, 1);  % 假设 Baselines 的第一列是 x
    y = Baselines(:, 2:3);  % 假设 Baselines 的第二到第五列是 y 的四个不同列%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 定义可选模型和高阶多项式模型的最大阶数
    numCols = size(y, 2);

    % 初始化存储最优模型的 cell 数组
    bestModel = cell(1, numCols);

    % 对每一列的 y 进行拟合，并选择最优模型
    for i = 1:numCols
        % 拟合数据
        % 使用内置的 gauss3 模型
        ft = fittype('gauss3');
        % startPoints = [0.6, 100, 40]; % 初始猜测值
        % startPoints = [0.6, 60, 20, 0.5, 140, 20]; % 初始猜测值 gauss2
        startPoints = [0.6, 60, 20, -0.4, 100, 20, 0.5, 140, 20]; % 初始猜测值 gauss3

        opts = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', startPoints, 'MaxIter', 1000, 'TolFun', 1e-9); % 调整函数容差

        fitResult = fit(x, y(:, i), ft, opts);
        % fitResult = fit(x, y(:, i), ft, opts, 'Start', startPoints);
        bestModel{i} = fitResult;

        % 显示最优模型的拟合参数
        fprintf('No %d column optimal model is: %s\n', i, formula(bestModel{i}));
    end
    %% %%%%%%%%%%%%%

    % 显示拟合参数
    disp('Fitted parameters:');
    disp(fitResult);

end
