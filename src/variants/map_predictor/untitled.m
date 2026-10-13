% 假设有 19 个肌肉的 EMG 数据，每个数据包含 101 个时间点
numMuscles = 19;
gaitCycle = linspace(0, 100, 101);

% 生成随机数据作为示例
experimentalEMG = rand(numMuscles, 101);
ssmMap = rand(numMuscles, 101);
sgmMap = rand(numMuscles, 101);

% 设置子图布局，5 列（对应不同百分比），每列 1 行
numColumns = 5;
numRows = ceil(numMuscles / numColumns);

% 创建图形窗口并设置尺寸
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

for i = 1:numMuscles
    % 创建子图
    subplot(numRows, numColumns, i);
    hold on;

    % 绘制实验数据
    plot(gaitCycle, experimentalEMG(i, :), '-k');

    % 绘制模型数据
    plot(gaitCycle, ssmMap(i, :), '--r');
    plot(gaitCycle, sgmMap(i, :), '--g');

    % 设置子图标题和标签
    if mod(i, numColumns) == 1
        ylabel('EMG Signal');
    end
    if i > (numRows - 1) * numColumns
        xlabel('Gait Cycle (%)');
    end
    title(sprintf('Muscle %d', i));
    hold off;
end

% 添加图例
legend({'Experimental EMG', 'SSM MAP', 'SGM MAP'}, 'Location', 'northeastoutside');
