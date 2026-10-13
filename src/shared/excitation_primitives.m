function [a, b, c, d, e, f, g, h] = excitation_primitives(Baselines)

    % 生成一些示例数据
    x = (1:200)';
    y = Baselines; % 加入一些噪声

    % 定义指数模型
    gaussModel = fittype('a*exp(-((x - b)^2) / (2*c^2)) + d + e*exp(-((x - f)^2) / (2*g^2)) + h', 'independent', 'x', 'coefficients', {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'});

    % 初始化存储拟合结果的结构数组
    fitResults = cell(1, 4);
    gofResults = cell(1, 4);
    fittedModels = cell(1, 4);

    % 进行拟合
    for i = 1:4

        [fitResults{i}, gofResults{i}] = fit(x, y(:, i), gaussModel);

        % extract the fitting results
        a(i) = fitResults{i}.a;
        b(i) = fitResults{i}.b;
        c(i) = fitResults{i}.c;
        d(i) = fitResults{i}.d;
        e(i) = fitResults{i}.e;
        f(i) = fitResults{i}.f;
        g(i) = fitResults{i}.g;
        h(i) = fitResults{i}.h;

        fittedModels{i} = @(x) a(i) * exp(-((x - b(i)).^2) / (2 * c(i).^2)) + d(i) + e(i) * exp(-((x - f(i)).^2) / (2 * g(i).^2)) + h(i);

    end

    x_new = (1:200)';  % 新的 x 值
    y_new = zeros(length(x_new), 4);
    for i = 1:4
        y_new(:, i) = fittedModels{i}(x_new);
    end

    figure;
    for i = 1:4
        subplot(2, 2, i);
        plot(x, y(:, i), 'o');  % 原始数据
        hold on;
        plot(x_new, y_new(:, i), '-r');  % 拟合曲线
        title(sprintf('第 %d 列的双指数函数拟合', i));
        xlabel('x');
        ylabel('y');
        legend('数据', '拟合曲线');
        hold off;
    end

end
