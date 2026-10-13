% LUMPED_MODEL_CHECK_WITH_COEFFS_THESIS Compare model predictions with the selected experimental reference signals.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

parentFolder = exo_data_root();
childFolder_1 = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_45\';
childFolder_2 = '\240528Pilot\CalibrationData\dynamics_stoop_15kg_45\';
[Reference_1, Scale_1] = Dynamic_kine_reference(childFolder_1, 1, 300);
[Reference_2, Scale_2] = Dynamic_kine_reference(childFolder_2, 154, 200);

start_stamp = 9033; % 1062 9033(5kg) 17033 8101(15kg) 9583(10kg)
end_stamp = 12248; % 4177、12248(5kg) 19809 10642(15kg) 12346(10kg)
start_stamp_2 = 17033; % 1062 9033 17033
end_stamp_2 = 19809; % 4177、12248 19809

%% \240730Exp\CalibrationData\dynamics_stoop_5kg_45_3\246

Torque_real_1 = Reference_1(:, 1);
Smooth_activation_1 = Reference_1(:, 2);
theta_real_1 = Reference_1(:, 3);
Target_LMT_1 = Reference_1(:, 4);

Torque_real_2 = Reference_2(:, 1);
Smooth_activation_2 = Reference_2(:, 2);
theta_real_2 = Reference_2(:, 3);
Target_LMT_2 = Reference_2(:, 4);

%% %%%%%%%%%%%%%%  edit the activations, theta, and torque to fit the dynamic mode-2 %%%%%%%%%%%

found = false;  % 用于标记是否已经找到了第一次满足条件的i
for i = 1:length(theta_real_1)
    if i <= start_stamp || i >= end_stamp
        if i <= start_stamp_2 || i >= end_stamp_2
            Torque_real_1(i) = 30;
            Smooth_activation_1(i) = 0;
            theta_real_1(i) = 0.5;
        end
    end
end

Torque_real_1 = smoothdata(Torque_real_1, "gaussian", 300);
Smooth_activation_1 = smoothdata(Smooth_activation_1, "gaussian", 300);
% Torque_real_1(10336:11177) = smoothdata(Torque_real_1(10336:11177), "gaussian", 1200);
% Smooth_activation_1(10336:11177) = smoothdata(Smooth_activation_1(10336:11177), "gaussian", 1200);
Torque_real_2 = smoothdata(Torque_real_2, "gaussian", 300);
Smooth_activation_2 = smoothdata(Smooth_activation_2, "gaussian", 300);
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

m1 = 5;

m2 = 15;

coeffs = [0.0400474589466669    0.700533596482596   0.175512162847502   10550.2360104555    45.8542898388227    0.286467647976371   0.798896309352907   1.04710175828713    1.89298082310009    54.8444188936951    0.700000000000000];

Zero_activation = zeros(length(Smooth_activation_2), 1)';
One_activation = ones(length(Smooth_activation_1), 1)';

% for i = 1:length(Smooth_activation_1)
%     if i < 15000 && i > 5000
%         One_activation(i) = 0.5;
%     end
% end
Model_torque = Torque_function_constant(m1, 0, theta_real_1, Torque_real_1, Smooth_activation_1, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11));

Model_torque_2 = Torque_function_constant(m1, 50, theta_real_1, Torque_real_1, Smooth_activation_1, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11));

Model_torque_3 = Torque_function_constant(m1, 0, theta_real_2, Torque_real_2, Zero_activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10), coeffs(11));

LMT = LMT_function(Smooth_activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the joint torque obtain by one rebuild muscle

figure;

plot(t_interp, Torque_real);

hold on;

plot(t_interp, Model_torque);

ylabel('L5S1_torque');

legend('Real-torque', 'Model-torque');
title('Real-torque, Model-torque');

figure;
plot(t_interp, activation);
hold on;

title('Average activation, EMGfilt-IL-r');

figure;
plot(t_interp, LMT);
title('LMT');
figure;
plot(t_interp, BE);
title('MA');

%% plot the model_LMT and target_LMT
figure;

plot(t_interp, Target_LMT);

hold on;

plot(t_interp, LMT);

ylabel('Muscle tendon length');

legend('Target-LMT', 'Model-LMT');
title('Target-LMT, Model-LMT');
