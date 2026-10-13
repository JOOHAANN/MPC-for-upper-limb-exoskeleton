%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% test_MAP.m
%%%
%%% Combine every code together and plot the MAP predictor
%%%
%%% Input: pv_coefficients.mat and pw_coefficients.mat
%%% Output: figures, MEP
%%%
%%% Youhan Li (2025)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%
% baseline_file = baseline_file(1:200);
%
% filename = 'baseline_data.mat';
%
% % save the data in the baseline_file
% save(filename, 'baseline_file');
load('baseline_data.mat');
load('EMG_matrix.mat');

load('Baselines.mat');

speed = [25, 35, 45, 55, 65];
weight = [5, 7.5, 10, 12.5, 15];

%% %%%%%%%%%%%
baseline_mode = 'SM';
muscle_titles = {'IL_r', 'IL_l', 'LTPT_r', 'LTPT_l', 'LTPL_r', 'LTPL_l'};

x = (1:200)';
y = Baselines; % 加入一些噪声

% 创建 Baselines 矩阵
Baselines = [x, y];

XP = getXP();

load('Primitives_fit_x.mat');
load('Primitives_fit_y.mat');

figure;
for i = 1:4
    subplot(2, 2, i);
    plot(x, y(:, i), 'o');  % raw data
    hold on;
    plot(Primitives_fit_x, Primitives_fit_y(:, i), '-r');  % fitting curve
    title(sprintf('number %d column of optimal combination', i));
    xlabel('x');
    ylabel('y');
    legend('data', 'fitting');
    % compute the correlated coeff
    R = corrcoef(Primitives_fit_y(:, i), y(:, i));
    r = R(1, 2);
    % transform the variable to the Str
    titleStr = sprintf('%s (r = %.2f)', r);

    % setting the sub title
    title(titleStr);
    hold off;
end

% j = 1;
%
% for ic = 1:5
%     for jc = 1:5
%
% [weights, XP, MEP] = MAP_model(speed(jc), weight(ic), baseline_mode, 'baseline_data.mat');
%
% figure
%
% k = 1:200;
%
% for i = 1:6
%     subplot(3, 2, i);
%     plot(k, MEP(:, i), 'b');
%     hold on;
%     plot(k, EMG_matrix(j).data(:, i), '--', 'Color', 'r');
%     legend('MEP', 'original data');
%
%     % compute the correlated coeff
%     R = corrcoef(MEP(:, i), EMG_matrix(j).data(:, i));
%     correlation(6 * (j - 1) + i) = R(1, 2);
%
%     errors = EMG_matrix(j).data(:, i) - MEP(:, i);
%     squared_errors = errors .^ 2;
%     mean_squared_error = mean(squared_errors);
%     rmse(6 * (j - 1) + i) = sqrt(mean_squared_error);
%
%     % transform the variable to the Str
%     titleStr = sprintf('%s (r = %.2f)', muscle_titles{i}, correlation(j));
%
%     % setting the sub title
%     title(titleStr);
% end
%
% titleStr = sprintf('dynamic %d kg %d', weight(ic), speed(jc));
% sgtitle(titleStr);
%
% j = j + 1;
%
%     end
% end

% 假设已定义:
%   speed  = [25,35,45,55,65];        % 5种速度
%   weight = [5,7,9,11,15];           % 5种重量
%   muscle_titles = {'IL_l','IL_r','LTPT_l','LTPT_r','LTPL_l','LTPL_r'};
%   (或其他类似的 6 个肌肉名，其中包含 _l 和 _r）

% 先筛选出只包含 "_r" 的肌肉索引
keepMask = contains(muscle_titles, '_r');
% 若你只想保留 IL_r, LTPT_r, LTPL_r，这里 keepMask 就是 [false true false true false true] 这样的逻辑向量
muscleIdx = find(keepMask);
% 这样 muscleIdx 只包含“带 _r”的那几个下标（通常是3个）

nMuscles = numel(muscleIdx); % 应该是 3

% 其余跟原先相同
% 还是要计算 correlation / rmse 的数组(如果需要)
correlation = zeros(5 * 5 * 3, 1); % 或者干脆只存 5*5*nMuscles 也行
rmse        = zeros(5 * 5 * 3, 1);

% 大循环 j 记录 EMG_matrix(j)
j = 1;

% -- 新建大图（只建一次），不再在循环里建 figure --
figure('Name', 'All-in-One', 'Color', 'w', 'DefaultAxesFontSize', 15, 'DefaultTextFontSize', 15);

% 我们要 15×5 个子图 (因为 5组重量×3肌肉=15 行，5组速度=5 列)
nRows = 5 * nMuscles;  % 5*3=15
nCols = 5;            % 5  (对应 5 种速度)

% 横坐标映射到 0~100
k = linspace(0, 100, 200);

for ic = 1:5  % weight
    for jc = 1:5  % speed

        % 调用你的模型函数
        [~, ~, MEP] = MAP_model(speed(jc), weight(ic), ...
                                baseline_mode, 'baseline_data.mat');

        % 遍历所有肌肉
        for m = 1:nMuscles
            % 计算“全局”行号
            rowIndex = nMuscles * (ic - 1) + m;

            % subplot 的索引
            idx = (rowIndex - 1) * nCols + jc;

            % 创建子图
            subplot(nRows, nCols, idx);
            hold on;

            % (1) 蓝色线 = 模型输出
            plot(k, MEP(:, m), 'b', 'LineWidth', 1.7);
            % (2) 红色线 = 原始数据
            plot(k, EMG_matrix(j).data(:, m), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.0);

            % 计算 correlation 和 RMSE
            Rmat = corrcoef(MEP(:, m), EMG_matrix(j).data(:, m));
            correlation(3 * (j - 1) + m) = Rmat(1, 2);

            errors = EMG_matrix(j).data(:, m) - MEP(:, m);
            rmse(3 * (j - 1) + m) = sqrt(mean(errors.^2));

            % 设置 X 轴
            xlim([0, 100]);
            if rowIndex == nRows
                set(gca, 'XTick', [0, 20, 40, 60, 80, 100]);
                set(gca, 'XTickLabel', {'0', '20', '40', '60', '80', '100'});
                xlabel('movement cycle (%)');
            else
                set(gca, 'XTick', []);
                set(gca, 'XTickLabel', '');
            end

            % 标注肌肉名称
            if jc == 1
                ylabel(muscle_titles{m});
                set(gca, 'YTick', []);
            end

            % 在第一行显示速度标题
            if rowIndex <= nMuscles && m == 1
                title(sprintf('%d%%', speed(jc)), 'FontWeight', 'bold');
            end

            % 在最左列显示重量
            if jc == 1 && m == ceil(nMuscles / 2)
                ylabel(sprintf('%.1f kg\n%s', weight(ic), muscle_titles{m}), 'FontSize', 12);
            end
        end

        % j 自增
        j = j + 1;

    end
end

% ----- 在大图右上角做一个总的 legend -----
% 先画“假线”只用于 legend
dummyBlue = plot(nan, nan, 'b', 'LineWidth', 1);
dummyRed  = plot(nan, nan, 'r', 'LineWidth', 1);
legend([dummyBlue, dummyRed], ...
       {'Model output', 'Reference data'}, ...
       'Location', 'northeastoutside');
%% plot IL_r, LTPT_r

% 获取包含 IL_r 和 LTPT_r 的肌肉索引
group1Mask = contains(muscle_titles, {'IL_r', 'LTPT_r'});
group1Idx = find(group1Mask);
nMuscles1 = numel(group1Idx);

% 创建第一个大图
figure('Name', 'Group1: IL_r & LTPT_r', 'Color', 'w', 'DefaultAxesFontSize', 18, 'DefaultTextFontSize', 18);

nRows1 = 5 * nMuscles1; % 行数
nCols = 5; % 列数 (速度)

j = 1; % 重置大循环计数器

for ic = 1:5  % weight
    for jc = 1:5  % speed
        [~, ~, MEP] = MAP_model(speed(jc), weight(ic), baseline_mode, 'baseline_data.mat');
        for m = 1:nMuscles1
            i = group1Idx(m);
            rowIndex = nMuscles1 * (ic - 1) + m;
            colIndex = jc;
            idx = (rowIndex - 1) * nCols + colIndex;
            subplot(nRows1, nCols, idx);
            hold on;
            plot(k, MEP(:, i), 'b', 'LineWidth', 1.7);
            plot(k, EMG_matrix(j).data(:, i), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.0);
            % Rmat = corrcoef(MEP(:, i), EMG_matrix(j).data(:, i));
            % correlation(6*(j-1) + i) = Rmat(1,2);
            % errors = EMG_matrix(j).data(:, i) - MEP(:, i);
            % rmse(6*(j-1) + i) = sqrt(mean(errors.^2));
            xlim([0, 100]);
            if rowIndex == nRows1
                set(gca, 'XTick', [0, 20, 40, 60, 80, 100]);
                xlabel('movement cycle (%)');
            else
                set(gca, 'XTick', []);
            end
            if colIndex == 1
                ylabel(muscle_titles{i});
                set(gca, 'YTick', []);
            end
            if (rowIndex <= nMuscles1) && (m == 1)
                title(sprintf('%d%', speed(jc)), 'FontWeight', 'bold');
            end
            if colIndex == 1 && m == ceil(nMuscles1 / 2)
                ylabel(sprintf('%.1f kg\n%s', weight(ic), muscle_titles{i}), 'FontSize', 18);
            end
        end
        j = j + 1;
    end
end

% 添加图例
dummyBlue = plot(nan, nan, 'b', 'LineWidth', 1);
dummyRed  = plot(nan, nan, 'r', 'LineWidth', 1);
legend([dummyBlue, dummyRed], {'Model output', 'Reference data'}, 'Location', 'northeastoutside');

%% plot LTPL_r

% 获取包含 LTPL_r 的肌肉索引
group2Mask = contains(muscle_titles, 'LTPL_r');
group2Idx = find(group2Mask);
nMuscles2 = numel(group2Idx);

% 创建第二个大图
figure('Name', 'Group2: LTPL_r', 'Color', 'w', 'DefaultAxesFontSize', 18, 'DefaultTextFontSize', 18);

nRows2 = 5 * nMuscles2; % 行数
nCols = 5; % 列数 (速度)

j = 1; % 重置大循环计数器

for ic = 1:5  % weight
    for jc = 1:5  % speed
        [~, ~, MEP] = MAP_model(speed(jc), weight(ic), baseline_mode, 'baseline_data.mat');
        for m = 1:nMuscles2
            i = group2Idx(m);
            rowIndex = nMuscles2 * (ic - 1) + m;
            colIndex = jc;
            idx = (rowIndex - 1) * nCols + colIndex;
            subplot(nRows2, nCols, idx);
            hold on;
            plot(k, MEP(:, i), 'b', 'LineWidth', 1.7);
            plot(k, EMG_matrix(j).data(:, i), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.0);
            % Rmat = corrcoef(MEP(:, i), EMG_matrix(j).data(:, i));
            % correlation(6*(j-1) + i) = Rmat(1,2);
            % errors = EMG_matrix(j).data(:, i) - MEP(:, i);
            % rmse(6*(j-1) + i) = sqrt(mean(errors.^2));
            xlim([0, 100]);
            if rowIndex == nRows2
                set(gca, 'XTick', [0, 20, 40, 60, 80, 100]);
                xlabel('movement cycle (%)');
            else
                set(gca, 'XTick', []);
            end
            if colIndex == 1
                ylabel(muscle_titles{i});
                set(gca, 'YTick', []);
            end
            if (rowIndex <= nMuscles2) && (m == 1)
                title(sprintf('%d%', speed(jc)), 'FontWeight', 'bold');
            end
            if colIndex == 1 && m == ceil(nMuscles2 / 2)
                ylabel(sprintf('%.1f kg\n%s', weight(ic), muscle_titles{i}), 'FontSize', 18);
            end
        end
        j = j + 1;
    end
end

% 添加图例
dummyBlue = plot(nan, nan, 'b', 'LineWidth', 1);
dummyRed  = plot(nan, nan, 'r', 'LineWidth', 1);
legend([dummyBlue, dummyRed], {'Model output', 'Reference data'}, 'Location', 'northeastoutside');

%% plot the correlation and RMSE

% correlation
mean_correlation = mean(correlation);
std_correlation = std(correlation);

figure('DefaultAxesFontSize', 27, 'DefaultTextFontSize', 27);
subplot(2, 2, 1);
histogram(correlation, 'FaceColor', 'b');
hold on;

hMean = xline(mean_correlation, '--r', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean_correlation));
hMean.FontSize = 22;

hMeanPlusSD = xline(mean_correlation + std_correlation, '--g', 'LineWidth', 2, 'Label', sprintf('Mean + SD = %.2f', mean_correlation + std_correlation));
hMeanPlusSD.FontSize = 22;

hMeanMinusSD = xline(mean_correlation - std_correlation, '--g', 'LineWidth', 2, 'Label', sprintf('Mean - SD = %.2f', mean_correlation - std_correlation));
hMeanMinusSD.FontSize = 22;

xlabel('Correlation');
ylabel('Frequency');
title('(a)');
legend('Correlation', 'Mean', 'Mean ± SD', 'Location', 'northwest');
set(legend, 'FontSize', 22);
% text(70, 2, sprintf('Mean = %.2f', mean_correlation), 'FontSize', 16, 'Color', 'r');
% text(70, 1.8, sprintf('SD = %.2f', std_correlation), 'FontSize', 16, 'Color', 'g');

% RMSE
mean_rmse = mean(rmse);
std_rmse = std(rmse);

subplot(2, 2, 2);
histogram(rmse, 'FaceColor', 'b');
hold on;

hMean = xline(mean_rmse, '--r', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean_rmse));
hMean.FontSize = 22; % 设置具体字体大小

hMeanPlusSD = xline(mean_rmse + std_rmse, '--g', 'LineWidth', 2, 'Label', sprintf('Mean + SD = %.2f', mean_rmse + std_rmse));
hMeanPlusSD.FontSize = 22;

hMeanMinusSD = xline(mean_rmse - std_rmse, '--g', 'LineWidth', 2, 'Label', sprintf('Mean - SD = %.2f', mean_rmse - std_rmse));
hMeanMinusSD.FontSize = 22;

xlabel('RMSE');
ylabel('Frequency');
title('(b)');
legend('RMSE', 'Mean', 'Mean ± SD');
set(legend, 'FontSize', 22);
% text(70, 2, sprintf('Mean = %.2f', mean_rmse), 'FontSize', 16, 'Color', 'r');
% text(70, 1.8, sprintf('SD = %.2f', std_rmse), 'FontSize', 16, 'Color', 'g');

%%%%%%%%%%%%%%%%%%%% Plot the experimental EMG data%%%%%%%%%%%%%%%%%%

% 创建第一个大图
figure('Name', 'Group1: IL_r & LTPT_r', 'Color', 'w', 'DefaultAxesFontSize', 18, 'DefaultTextFontSize', 18);

nRows1 = 5 * nMuscles1; % 行数
nCols = 5; % 列数 (速度)

j = 1; % 重置大循环计数器

for ic = 1:5  % weight
    for jc = 1:5  % speed
        [~, ~, MEP] = MAP_model(speed(jc), weight(ic), baseline_mode, 'baseline_data.mat');
        for m = 1:nMuscles1
            i = group1Idx(m);
            rowIndex = nMuscles1 * (ic - 1) + m;
            colIndex = jc;
            idx = (rowIndex - 1) * nCols + colIndex;
            subplot(nRows1, nCols, idx);
            hold on;
            plot(k, EMG_matrix(j).data(:, i), 'Color', 'b', 'LineWidth', 1.8);
            % Rmat = corrcoef(MEP(:, i), EMG_matrix(j).data(:, i));
            % correlation(6*(j-1) + i) = Rmat(1,2);
            % errors = EMG_matrix(j).data(:, i) - MEP(:, i);
            % rmse(6*(j-1) + i) = sqrt(mean(errors.^2));
            xlim([0, 100]);
            if rowIndex == nRows1
                set(gca, 'XTick', [0, 20, 40, 60, 80, 100]);
                xlabel('movement cycle (%)');
            else
                set(gca, 'XTick', []);
            end
            if colIndex == 1
                ylabel(muscle_titles{i});
                set(gca, 'YTick', []);
            end
            if (rowIndex <= nMuscles1) && (m == 1)
                title(sprintf('%d%', speed(jc)), 'FontWeight', 'bold');
            end
            if colIndex == 1 && m == ceil(nMuscles1 / 2)
                ylabel(sprintf('%.1f kg\n%s', weight(ic), muscle_titles{i}), 'FontSize', 18);
            end
        end
        j = j + 1;
    end
end

% 添加图例
dummyBlue = plot(nan, nan, 'b', 'LineWidth', 1);
legend([dummyBlue, dummyRed], {'Experimental EMG data'}, 'Location', 'northeastoutside');

%% plot LTPL_r

% 获取包含 LTPL_r 的肌肉索引
group2Mask = contains(muscle_titles, 'LTPL_r');
group2Idx = find(group2Mask);
nMuscles2 = numel(group2Idx);

% 创建第二个大图
figure('Name', 'Group2: LTPL_r', 'Color', 'w', 'DefaultAxesFontSize', 18, 'DefaultTextFontSize', 18);

nRows2 = 5 * nMuscles2; % 行数
nCols = 5; % 列数 (速度)

j = 1; % 重置大循环计数器

for ic = 1:5  % weight
    for jc = 1:5  % speed
        [~, ~, MEP] = MAP_model(speed(jc), weight(ic), baseline_mode, 'baseline_data.mat');
        for m = 1:nMuscles2
            i = group2Idx(m);
            rowIndex = nMuscles2 * (ic - 1) + m;
            colIndex = jc;
            idx = (rowIndex - 1) * nCols + colIndex;
            subplot(nRows2, nCols, idx);
            hold on;
            plot(k, EMG_matrix(j).data(:, i), 'Color', 'b', 'LineWidth', 1.8);
            % Rmat = corrcoef(MEP(:, i), EMG_matrix(j).data(:, i));
            % correlation(6*(j-1) + i) = Rmat(1,2);
            % errors = EMG_matrix(j).data(:, i) - MEP(:, i);
            % rmse(6*(j-1) + i) = sqrt(mean(errors.^2));
            xlim([0, 100]);
            if rowIndex == nRows2
                set(gca, 'XTick', [0, 20, 40, 60, 80, 100]);
                xlabel('movement cycle (%)');
            else
                set(gca, 'XTick', []);
            end
            if colIndex == 1
                ylabel(muscle_titles{i});
                set(gca, 'YTick', []);
            end
            if (rowIndex <= nMuscles2) && (m == 1)
                title(sprintf('%d%', speed(jc)), 'FontWeight', 'bold');
            end
            if colIndex == 1 && m == ceil(nMuscles2 / 2)
                ylabel(sprintf('%.1f kg\n%s', weight(ic), muscle_titles{i}), 'FontSize', 18);
            end
        end
        j = j + 1;
    end
end

% 添加图例
dummyBlue = plot(nan, nan, 'b', 'LineWidth', 1);
legend([dummyBlue, dummyRed], {'Experimental EMG data'}, 'Location', 'northeastoutside');

%%%%%%%%%%%%%%%%%%%% Plot the experimental EMG data%%%%%%%%%%%%%%%%%%
