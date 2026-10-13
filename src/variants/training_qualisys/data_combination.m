% DATA_COMBINATION Assemble selected movement segments into the calibration data matrix.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

%%%%%%% step2%%%%%%%% select the segment
parentFolder = exo_data_root();
childFolder_1 = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_45\';
childFolder_2 = '\240528Pilot\CalibrationData\dynamics_stoop_10kg_35\';
childFolder_3 = '\240528Pilot\CalibrationData\dynamics_stoop_15kg_45\';
%%%%%%% step2%%%%%%%% select the segment

%%%%%%% step2%%%%%%%% select the segment

[Reference_1, Scale_1] = Dynamic_kine_reference(childFolder_1, 1, 300);
[Reference_2, Scale_2] = Dynamic_kine_reference(childFolder_2, 1, 350);
[Reference_3, Scale_3] = Dynamic_kine_reference(childFolder_3, 50, 300);

%%%%%%% step2%%%%%%%% select the segment

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

%% %%%%%%%%%%%%%%  edit the activations, theta, and torque to fit the dynamic mode-2 %%%%%%%%%%%

%%%%%%%%%%%% sewing the model%%%%%%%%%%%%

Torque_sew = [Torque_real_1; Torque_real_2; Torque_real_3];
Smooth_activation_sew = [Smooth_activation_1; Smooth_activation_2; Smooth_activation_3];
theta_real_sew = [theta_real_1; theta_real_2; theta_real_3];
Smooth_activation_front_sew = [Smooth_activation_front_1; Smooth_activation_front_2; Smooth_activation_front_3];

%%%%%%% step3%%%%%%%% select the segment

start_stamp_1 = 7757; % 1062 9033(5kg) 17033 8101(15kg) 9583(10kg)
end_stamp_1 = 13853; % 4177、12248(5kg) 19809 10642(15kg) 12346(10kg)
start_stamp_2 = 20956; % 1062 9033 17033
end_stamp_2 = 27639; % 4177、12248 19809
start_stamp_3 = 36908; % 1062 9033(5kg) 17033 8101(15kg) 9583(10kg)
end_stamp_3 = 41734; % 4177、12248(5kg) 19809 10642(15kg) 12346(10kg)

%%%%%%% step3%%%%%%%% select the segment

%%%%%%%%%%%% sewing the model%%%%%%%%%%%%

found = false;
for i = 1:length(Torque_sew)
    if i <= start_stamp_1 || i >= end_stamp_1 && i <= start_stamp_2 || i >= end_stamp_2 && i <= start_stamp_3 || i >= end_stamp_3
        Torque_sew(i) = 20;
        Smooth_activation_sew(i) = 0;
        Smooth_activation_front_sew(i) = 0;
        theta_real_sew(i) = 1.5;
    end
end

theta_real_sew = smoothdata(theta_real_sew, "gaussian", 300);
Torque_sew = smoothdata(Torque_sew, "gaussian", 300);
Smooth_activation_sew = smoothdata(Smooth_activation_sew, "gaussian", 300);
Smooth_activation_front_sew = lowpassfilter(Smooth_activation_front_sew, 1000, 1 / 1.2);

Combination(:, 1) = Torque_sew;
Combination(:, 2) = Smooth_activation_sew;
Combination(:, 3) = theta_real_sew;
Combination(:, 4) = Smooth_activation_front_sew;
save('Combination.mat', 'Combination');
