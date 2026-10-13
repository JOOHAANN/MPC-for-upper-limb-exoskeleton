% DATA_COMBINATION Assemble selected movement segments into the calibration data matrix.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

parentFolder = exo_data_root();
childFolder_1 = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_45\';
childFolder_2 = '\240528Pilot\CalibrationData\dynamics_stoop_10kg_35\';
childFolder_3 = '\240528Pilot\CalibrationData\dynamics_stoop_15kg_55\';
childFolder_4 = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_25\';
childFolder_5 = '\240528Pilot\CalibrationData\dynamics_stoop_10kg_25\';
childFolder_6 = '\240528Pilot\CalibrationData\dynamics_stoop_15kg_25\';
childFolder_7 = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_65\';
childFolder_8 = '\240528Pilot\CalibrationData\dynamics_stoop_10kg_65\';
childFolder_9 = '\240528Pilot\CalibrationData\dynamics_stoop_15kg_65\';
[Reference_1, Scale_1] = Dynamic_kine_reference(childFolder_1, 1, 300);
[Reference_2, Scale_2] = Dynamic_kine_reference(childFolder_2, 1, 350);
[Reference_3, Scale_3] = Dynamic_kine_reference(childFolder_3, 50, 300);
[Reference_4, Scale_4] = Dynamic_kine_reference(childFolder_4, 1, 300);
[Reference_5, Scale_5] = Dynamic_kine_reference(childFolder_5, 1, 300);
[Reference_6, Scale_6] = Dynamic_kine_reference(childFolder_6, 1, 300);
[Reference_7, Scale_7] = Dynamic_kine_reference(childFolder_7, 1, 300);
[Reference_8, Scale_8] = Dynamic_kine_reference(childFolder_8, 1, 300);
[Reference_9, Scale_9] = Dynamic_kine_reference(childFolder_9, 1, 300);

start_stamp_4 = 17033; % 1062 9033 17033
end_stamp_4 = 19809; % 4177、12248 19809
start_stamp_5 = 7532; % 1062 9033(5kg) 17033 8101(15kg) 9583(10kg)
end_stamp_5 = 11203; % 4177、12248(5kg) 19809 10642(15kg) 12346(10kg)
start_stamp_6 = 17033; % 1062 9033 17033
end_stamp_6 = 19809; % 4177、12248 19809
start_stamp_7 = 7532; % 1062 9033(5kg) 17033 8101(15kg) 9583(10kg)
end_stamp_7 = 11203; % 4177、12248(5kg) 19809 10642(15kg) 12346(10kg)
start_stamp_8 = 17033; % 1062 9033 17033
end_stamp_8 = 19809; % 4177、12248 19809
start_stamp_9 = 17033; % 1062 9033 17033
end_stamp_9 = 19809; % 4177、12248 19809

%% \240730Exp\CalibrationData\dynamics_stoop_5kg_45_3\246

Torque_real_1 = Reference_1(:, 1);
Smooth_activation_1 = Reference_1(:, 2);
theta_real_1 = Reference_1(:, 3);
Target_LMT_1 = Reference_1(:, 4);
Smooth_activation_front_1 = Reference_1(:, 5);

Torque_real_2 = Reference_2(:, 1);
Smooth_activation_2 = Reference_2(:, 2);
theta_real_2 = Reference_2(:, 3);
Target_LMT_2 = Reference_2(:, 4);
Smooth_activation_front_2 = Reference_2(:, 5);

Torque_real_3 = Reference_3(:, 1);
Smooth_activation_3 = Reference_3(:, 2);
theta_real_3 = Reference_3(:, 3);
Target_LMT_3 = Reference_3(:, 4);
Smooth_activation_front_3 = Reference_3(:, 5);

Torque_real_4 = Reference_4(:, 1);
Smooth_activation_4 = Reference_4(:, 2);
theta_real_4 = Reference_4(:, 3);
Target_LMT_4 = Reference_4(:, 4);

%% %%%%%%%%%%%%%%  edit the activations, theta, and torque to fit the dynamic mode-2 %%%%%%%%%%%

%%%%%%%%%%%% sewing the model%%%%%%%%%%%%

Torque_sew = [Torque_real_1; Torque_real_2; Torque_real_3];
Smooth_activation_sew = [Smooth_activation_1; Smooth_activation_2; Smooth_activation_3];
theta_real_sew = [theta_real_1; theta_real_2; theta_real_3];
Target_LMT_sew = [Target_LMT_1; Target_LMT_2; Target_LMT_3];

Smooth_activation_front_sew = [Smooth_activation_front_1; Smooth_activation_front_2; Smooth_activation_front_3];

start_stamp_1 = 7775; % 1062 9033(5kg) 17033 8101(15kg) 9583(10kg)
end_stamp_1 = 13388; % 4177、12248(5kg) 19809 10642(15kg) 12346(10kg)
start_stamp_2 = 21480; % 1062 9033 17033
end_stamp_2 = 27486; % 4177、12248 19809
start_stamp_3 = 39126; % 1062 9033(5kg) 17033 8101(15kg) 9583(10kg)
end_stamp_3 = 42890; % 4177、12248(5kg) 19809 10642(15kg) 12346(10kg)

%%%%%%%%%%%% sewing the model%%%%%%%%%%%%

found = false;  % 用于标记是否已经找到了第一次满足条件的i
for i = 1:length(Torque_sew)
    if i <= start_stamp_1 || i >= end_stamp_1 && i <= start_stamp_2 || i >= end_stamp_2 && i <= start_stamp_3 || i >= end_stamp_3
        Torque_sew(i) = 20;
        Smooth_activation_sew(i) = 0;
        Smooth_activation_front_sew(i) = 0;
        theta_real_sew(i) = 1.5;
        Target_LMT_sew(i) = 0.222;
    end
end

theta_real_sew = smoothdata(theta_real_sew, "gaussian", 300);
Target_LMT_sew = smoothdata(Target_LMT_sew, "gaussian", 300);
Torque_sew = smoothdata(Torque_sew, "gaussian", 300);
Smooth_activation_sew = smoothdata(Smooth_activation_sew, "gaussian", 300);
Smooth_activation_front_sew = lowpassfilter(Smooth_activation_front_sew, 1000, 1 / 1.2);

% Torque_real_1(10336:11177) = smoothdata(Torque_real_1(10336:11177), "gaussian", 1200);
% Smooth_activation_1(10336:11177) = smoothdata(Smooth_activation_1(10336:11177), "gaussian", 1200);
Combination(:, 1) = Torque_sew;
Combination(:, 2) = Smooth_activation_sew;
Combination(:, 3) = theta_real_sew;
Combination(:, 4) = Target_LMT_sew;
save('Combination_raw.mat', 'Combination');
