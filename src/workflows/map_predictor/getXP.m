%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% getXP.mat
%%%
%%% Estimates the exitacion primitive (XP) for a specified speed and
%%% elevation. There are 2 baseline modes: 'GM', which is the generic mode
%%% and no file should be indicated and 'SM', which is the subject-specific
%%% mode and a file should be indicated in the next argument (baseline_file)
%%%
%%% Jose Gonzalez (2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XP] = getXP()

    load('Baselines.mat');

    % [a,b,c,d,e,f,g,h] = excitation_primitives(Baselines);

    %% %%%%%%%%%%%%%

    x = (1:200)';
    y = Baselines; % 加入一些噪声

    % 创建 Baselines 矩阵
    Baselines = [x, y];

    % 调用 auto_select_model 函数
    bestModels = excitation_primitives_double_gauss(Baselines);

    % 使用提取的最优模型计算新的 y 值
    Primitives_fit_x = (1:200)';  % 新的 x 值
    Primitives_fit_y = zeros(length(Primitives_fit_x), 4);
    for i = 1:4
        Primitives_fit_y(:, i) = feval(bestModels{i}, Primitives_fit_x);
    end

    filename_Primitives_fit_x = 'Primitives_fit_x.mat';
    filename_Primitives_fit_y = 'Primitives_fit_y.mat';

    save(filename_Primitives_fit_x, 'Primitives_fit_x');
    save(filename_Primitives_fit_y, 'Primitives_fit_y');

    % 绘制结果
    % figure;
    % for i = 1:4
    %     subplot(2, 2, i);
    %     plot(x, y(:, i), 'o');  % raw data
    %     hold on;
    %     plot(Primitives_fit_x, Primitives_fit_y(:, i), '-r');  % fitting curve
    %     title(sprintf('number %d column of optimal combination', i));
    %     xlabel('x');
    %     ylabel('y');
    %     legend('data', 'fitting');
    %     hold off;
    % end

    % the fitted excitation primitives curves

    XP = Primitives_fit_y;

end
