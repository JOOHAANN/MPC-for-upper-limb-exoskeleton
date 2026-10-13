% LUMPED_MODEL_CHECK_WITH_COEFFS_ICNR Compare model predictions with the selected experimental reference signals.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

parentFolder = exo_data_root();
childFolder_1 = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_45\';
childFolder_2 = '\240528Pilot\CalibrationData\dynamics_stoop_10kg_45\';
childFolder_3 = '\240528Pilot\CalibrationData\dynamics_stoop_15kg_45\';
Reference_1 = Dynamic_kine_reference(childFolder_1, 163, 253); %% 163, 253
Reference_2 = Dynamic_kine_reference(childFolder_2, 183, 273);
Reference_3 = Dynamic_kine_reference(childFolder_3, 137, 227);

Torque_real_1 = smoothdata(Reference_1(:, 1), "gaussian", 5);
Smooth_activation_1 = Reference_1(:, 2);
theta_real_1 = Reference_1(:, 3);
Target_LMT_1 = Reference_1(:, 4);

Torque_real_2 = smoothdata(Reference_2(:, 1), "gaussian", 5);
Smooth_activation_2 = Reference_2(:, 2);
theta_real_2 = Reference_2(:, 3);
Target_LMT_2 = Reference_2(:, 4);

Torque_real_3 = smoothdata(Reference_3(:, 1), "gaussian", 5);
Smooth_activation_3 = Reference_3(:, 2);
theta_real_3 = Reference_3(:, 3);
Target_LMT_3 = Reference_3(:, 4);

Target_theta = cell(3, 1);
Target_theta{1} = theta_real_1;
Target_theta{2} = theta_real_2;
Target_theta{3} = theta_real_3;

Target_torque = cell(3, 1);
Target_torque{1} = Torque_real_1;
Target_torque{2} = Torque_real_2;
Target_torque{3} = Torque_real_3;

%% Exp240405
%% 240408Youhan
%% \240319Pilot\dynamic_stoop_5kg\
%% \stoop_5kg_3\
%% \240430Youhan\CalibrationData\pilot_noexo_stoop_5_2\
%% \240507Pilot\dynamic_stoop_5kg_3\
%% \240528Pilot\CalibrationData\dynamics_stoop_15kg_45\
%% \240514Pilot\CalibrationData\dynamic_stoop_5_0\
%% 240725Pilot\dynamics_stoop_5kg
%% 240730Exp\CalibrationData\dynamics_stoop_5kg_45_3
Ts = 0.001;

damping = 0.5;

alpha0 = 0.24086; % obtain from the paper
% 1924.6

m1 = 5;

m2 = 10;

m3 = 15;

coeffs = [0.0307584358406881    0.337948403010836   0.304276516510970   6911.26423096057    20.5156249020510    0.762060141207511   0.490398273940264   0.275995566721535   0.700000000000000];

Model_torque_1 = Torque_function_constant(m1, 0, theta_real_1, Smooth_activation_1, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9));

Model_torque_2 = Torque_function_constant(m2, 100, theta_real_1, Smooth_activation_1, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9));

Model_torque_3 = Torque_function_constant(m3, 0, theta_real_3, Smooth_activation_3, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9));

% 创建一个单元数组来存储不同尺寸的矩阵
Model_torque = cell(3, 1);

Model_torque{1} = Model_torque_1;
Model_torque{2} = Model_torque_2;
Model_torque{3} = Model_torque_3;
% [LMT, BE] = simplify_pulley_model(L, coeffs(1), Model_torque(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the joint torque obtain by one rebuild muscle

% figure;
%
% plot(t_interp - t_start, Torque_real);
%
% hold on;
%
% plot(t_interp - t_start, Model_torque(:,1));
%
% ylabel('L5S1_torque')
%
% legend('Real-torque', 'Model-torque');
% title('Real-torque, Model-torque');
%
% figure;
% plot(t_interp - t_start, activation);
% hold on;
% % plot(fmt.data(206:758,1) - fmt.data(206,1), EMGGG);
% % legend('Average activation', 'EMGfilt-IL-r');
% title('Average activation, EMGfilt-IL-r');
%
%
% %% plot the model_theta and target_theta
% figure;
%
% plot(t_interp - t_start, Target_theta);
%
% hold on;
%
% plot(t_interp - t_start, Model_theta);
%
% ylabel('Muscle tendon length')
%
% legend('Target_theta', 'Model_theta');
% title('Target_theta, Model_theta');

% %% plot the moment of arm
%
% figure;
% plot(t_interp - t_start, BE);
% title('MA');

% %% plot the all the outputs together
% % torque
% R = corrcoef(Model_torque(:,1), Torque_real);
% correlation(1) = R(1, 2);
% figure('DefaultAxesFontSize', 24, 'DefaultTextFontSize', 24);
% subplot(2,1,1);
% plot(t_interp - t_start, Torque_real, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3.5);
%
% hold on;
%
% plot(t_interp - t_start, Model_torque(:,1), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5);
% xlabel('Time (s)');
% ylabel('Torque(N.M)')
%
% Model_torque_max = max(Model_torque(:,1)); % rescale the profile
% ylim([0, Model_torque_max] * 1.3);
% legend('CEINMS-Torque', 'Model-Torque','FontSize', 20);
% title('(b)                   ');
% legend boxoff;
% % % MTF
% % subplot(2,1,2);
% % plot(t_interp - t_start, Model_MTF, 'LineWidth', 2);
% % xlabel('Time (s)');
% % ylabel('Force(N)')
% %
% % Model_MTF_max = max(Model_MTF); % rescale the profile
% % ylim([0, Model_MTF_max] * 1.1);
% %
% % legend('MTF','FontSize', 15);
% % title('(c) musculotendon force generate by the model');
%
% % LMT
%
% R = corrcoef(LMT, interp1(lmt.data(:,1), Target_LMT, t_interp));
% correlation(2) = R(1, 2);
% subplot(2,1,2);
%
%
% plot(t_interp - t_start, interp1(lmt.data(:,1), Target_LMT, t_interp), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3.5);
%
% hold on;
%
% plot(t_interp - t_start, LMT, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5);
%
% xlabel('Time (s)');
% ylabel('L^{MT}(m)')
% LMT_max = max(LMT) - min(Target_LMT); % rescale the profile
% LMT_min = min(Target_LMT); % rescale the profile
% ylim([LMT_min, LMT_min + LMT_max * 1.5]);
% legend('CEINMS-L^{MT}', 'Model-L^{MT}','FontSize', 20);
% title('(c)                      ');
% legend boxoff;
% % % moment arm
% %
% % subplot(2,2,4);
% % plot(t_interp - t_start, BE, 'LineWidth', 2);
% % xlabel('Time (s)');
% %
% % BE_max = max(BE) - min(BE); % rescale the profile
% % BE_min = min(BE); % rescale the profile
% % ylim([BE_min, BE_min + BE_max * 1.1]);
% %
% % legend('MA','FontSize', 15);
% % ylabel('MA(m)')
% % title('(d) Moment arm of the model');
%
% % 调整子图的位置和宽度，使其保持一致
% subplot_handles = findall(gcf, 'Type', 'Axes');
% for i = 1:length(subplot_handles)
%     subplot_handles(i).Position(3) = 0.3; % 设置子图的宽度
% end

%% plot the all the outputs together
% 创建一个 3x2 的子图网格
figure('Units', 'pixels', ...
       'DefaultAxesFontSize', 14, ...
       'DefaultTextFontSize', 14);

for i = 1:3
    % 获取当前的 torque 和 theta 数据

    MT = Model_torque{i};
    TT = Target_torque{i};
    TA = Target_theta{i};
    labels = ["(b)", "(c)", "(d)", "(e)", "(f)", "(g)"];
    % 绘制 torque 图
    R = corrcoef(MT(:, 1), TT);
    Torque_correlation(i) = R(1, 2);
    subplot(3, 2, 2 * i - 1);
    plot(TT, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3.5);
    hold on;
    plot(MT(:, 1), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5);
    xlabel('                                        Time (ms)');
    ylabel('Torque (N.M)');
    Model_torque_max = max(MT(:, 1)); % rescale the profile
    ylim([0, Model_torque_max] * 1.3);
    legend('CEINMS-Torque', 'Model-Torque', 'FontSize', 14);
    title(labels(2 * i - 1), 'FontSize', 22);
    legend boxoff;

    % 绘制 Theta 图
    R = corrcoef(MT(:, 2), TA);
    Theta_correlation(i) = R(1, 2);
    subplot(3, 2, 2 * i);
    plot(TA, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3.5);
    hold on;
    plot(MT(:, 2), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5);
    xlabel('                                        Time (ms)');
    ylabel('Theta (rad)');
    LMT_max = max(MT(:, 2)) - min(MT(:, 2)); % rescale the profile
    LMT_min = min(MT(:, 2)); % rescale the profile
    ylim([LMT_min, LMT_min + LMT_max * 1.5]);
    legend('10*CEINMS-theta^{MT}', 'Model-theta^{MT}', 'FontSize', 14);
    title(labels(2 * i), 'FontSize', 22);
    legend boxoff;
end

% 调整子图的位置和宽度
subplot_handles = findall(gcf, 'Type', 'Axes');
for i = 1:length(subplot_handles)
    % 设置每个子图的宽度和高度
    subplot_handles(i).Position(3) = 0.22; % 设置子图的宽度
    subplot_handles(i).Position(4) = 0.24; % 设置子图的高度
    if mod(i, 2) == 1
        subplot_handles(i).Position(1) = 0.4;
    end
end

% for i = 1:3
%
% % torque
% MT = Model_torque{i};
% TT = Target_torque{i};
% % R(i,:,:)  = corrcoef(MT(:,1), TT);
% % correlation(1) = R(:, 1, 2);
% figure('Units', 'pixels', ...
%        'Position', [100, 100, 600, 600], ...
%        'DefaultAxesFontSize', 24, ...
%        'DefaultTextFontSize', 24);
% subplot(3, 2, 2*i-1);
% plot(TT, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3.5);
%
% hold on;
%
% plot(MT(:,1), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5);
% xlabel('Time (s)');
% ylabel('Torque(N.M)')
%
% Model_torque_max = max(MT(:,1)); % rescale the profile
% ylim([0, Model_torque_max] * 1.3);
% legend('CEINMS-Torque', 'Model-Torque','FontSize', 20);
% title('(b)                   ');
% legend boxoff;
% % % MTF
% % subplot(2,1,2);
% % plot(t_interp - t_start, Model_MTF, 'LineWidth', 2);
% % xlabel('Time (s)');
% % ylabel('Force(N)')
% %
% % Model_MTF_max = max(Model_MTF); % rescale the profile
% % ylim([0, Model_MTF_max] * 1.1);
% %
% % legend('MTF','FontSize', 15);
% % title('(c) musculotendon force generate by the model');
%
% % LMT
%
% TA = Target_theta{i};
% % R(i,:,:) = corrcoef(MT(:,1), TA);
% % correlation(2) = R(i, 1, 2);
% subplot(3, 2, 2*i);
%
%
% plot(TA, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3.5);
%
% hold on;
%
% plot(MT(:,2), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5);
%
% xlabel('Time (s)');
% ylabel('L^{MT}(m)')
% LMT_max = max(MT(:,2)) - min(MT(:,2)); % rescale the profile
% LMT_min = min(MT(:,2)); % rescale the profile
% ylim([LMT_min, LMT_min + LMT_max * 1.5]);
% legend('CEINMS-L^{MT}', 'Model-L^{MT}','FontSize', 20);
% title('(c)                      ');
% legend boxoff;
%
% end
% % % moment arm
% %
% % subplot(2,2,4);
% % plot(t_interp - t_start, BE, 'LineWidth', 2);
% % xlabel('Time (s)');
% %
% % BE_max = max(BE) - min(BE); % rescale the profile
% % BE_min = min(BE); % rescale the profile
% % ylim([BE_min, BE_min + BE_max * 1.1]);
% %
% % legend('MA','FontSize', 15);
% % ylabel('MA(m)')
% % title('(d) Moment arm of the model');
%
% % 调整子图的位置和宽度，使其保持一致
% subplot_handles = findall(gcf, 'Type', 'Axes');
% for i = 1:length(subplot_handles)
%     subplot_handles(i).Position(3) = 0.3; % 设置子图的宽度
% end
