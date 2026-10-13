clc;
clear;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

parentFolder = exo_data_root();
childFolder = '\240528Pilot\CalibrationData\dynamics_stoop_5kg_45\';

%% Exp240405
%% 240408Youhan
%% \240319Pilot\dynamic_stoop_5kg\
%% \stoop_5kg_3\
%% \240430Youhan\CalibrationData\pilot_noexo_stoop_5_2\
%% \240507Pilot\dynamic_stoop_5kg_3\
%% \240528Pilot\CalibrationData\dynamics_stoop_15kg_45\
%% \240514Pilot\CalibrationData\dynamic_stoop_5_0\

EMGfilt = importdata(exo_data_file(parentFolder, childFolder, 'emgFilt.sto'));
ikData = importdata(exo_data_file(parentFolder, childFolder, 'ik.sto'));

lmt = importdata(exo_data_file(parentFolder, childFolder, 'lmt.sto'));
fmt = importdata(exo_data_file(parentFolder, childFolder, 'MuscleForces.sto'));
lm = importdata(exo_data_file(parentFolder, childFolder, 'FibreLengths.sto'));
lt = importdata(exo_data_file(parentFolder, childFolder, 'TendonLength.sto'));
vm = importdata(exo_data_file(parentFolder, childFolder, 'FibreVelocities.sto'));
ma = importdata(exo_data_file(parentFolder, childFolder, 'ma_L5_S1_Flex_Ext.sto'));
Tq = importdata(exo_data_file(parentFolder, childFolder, 'Torque.sto'));
Ac = importdata(exo_data_file(parentFolder, childFolder, 'Activations.sto'));

%% EMG Data

figure;
for iE = 2:length(EMGfilt.colheaders)
    subplot(4, 2, iE - 1);
    plot(EMGfilt.data(:, 1) - EMGfilt.data(1, 1), EMGfilt.data(:, iE));

    title([EMGfilt.colheaders{iE}]);
    ylabel('EMG (a.u)');
    ylim([0 1.1]);
end

%% L5S1 flexion-extension joint angles

figure;
plot(ikData.data(:, 1) - ikData.data(1, 1), ikData.data(:, 2) * (180 / pi));
% ylabel('L5S1 Joint Angle'); xlim([0,ikData.data(end,1)-ikData.data(1,1)])

%% Indexes identification

rect_abd_r_ind = find(contains(lmt.colheaders, 'rect_abd_r'));
rect_abd_l_ind = find(contains(lmt.colheaders, 'rect_abd_l'));

IL_L1_r_ind = find(contains(lmt.colheaders, 'IL_L1_r'));
IL_L2_r_ind = find(contains(lmt.colheaders, 'IL_L2_r'));
IL_L3_r_ind = find(contains(lmt.colheaders, 'IL_L3_r'));
IL_L4_r_ind = find(contains(lmt.colheaders, 'IL_L4_r'));

IL_L1_l_ind = find(contains(lmt.colheaders, 'IL_L1_l'));
IL_L2_l_ind = find(contains(lmt.colheaders, 'IL_L2_l'));
IL_L3_l_ind = find(contains(lmt.colheaders, 'IL_L3_l'));
IL_L4_l_ind = find(contains(lmt.colheaders, 'IL_L4_l'));

LTpT_T7_r_ind = find(contains(lmt.colheaders, 'LTpT_T7_r'));
LTpT_T8_r_ind = find(contains(lmt.colheaders, 'LTpT_T8_r'));
LTpT_T9_r_ind = find(contains(lmt.colheaders, 'LTpT_T9_r'));
LTpT_T10_r_ind = find(contains(lmt.colheaders, 'LTpT_T10_r'));
LTpT_T11_r_ind = find(contains(lmt.colheaders, 'LTpT_T11_r'));
LTpT_T12_r_ind = find(contains(lmt.colheaders, 'LTpT_T12_r'));

LTpT_R7_r_ind = find(contains(lmt.colheaders, 'LTpT_R7_r'));
LTpT_R8_r_ind = find(contains(lmt.colheaders, 'LTpT_R8_r'));
LTpT_R9_r_ind = find(contains(lmt.colheaders, 'LTpT_R9_r'));
LTpT_R10_r_ind = find(contains(lmt.colheaders, 'LTpT_R10_r'));
LTpT_R11_r_ind = find(contains(lmt.colheaders, 'LTpT_R11_r'));
LTpT_R12_r_ind = find(contains(lmt.colheaders, 'LTpT_R12_r'));

LTpT_T7_l_ind = find(contains(lmt.colheaders, 'LTpT_T7_l'));
LTpT_T8_l_ind = find(contains(lmt.colheaders, 'LTpT_T8_l'));
LTpT_T9_l_ind = find(contains(lmt.colheaders, 'LTpT_T9_l'));
LTpT_T10_l_ind = find(contains(lmt.colheaders, 'LTpT_T10_l'));
LTpT_T11_l_ind = find(contains(lmt.colheaders, 'LTpT_T11_l'));
LTpT_T12_l_ind = find(contains(lmt.colheaders, 'LTpT_T12_l'));

LTpT_R7_l_ind = find(contains(lmt.colheaders, 'LTpT_R7_l'));
LTpT_R8_l_ind = find(contains(lmt.colheaders, 'LTpT_R8_l'));
LTpT_R9_l_ind = find(contains(lmt.colheaders, 'LTpT_R9_l'));
LTpT_R10_l_ind = find(contains(lmt.colheaders, 'LTpT_R10_l'));
LTpT_R11_l_ind = find(contains(lmt.colheaders, 'LTpT_R11_l'));
LTpT_R12_l_ind = find(contains(lmt.colheaders, 'LTpT_R12_l'));

LTpL_L1_r_ind = find(contains(lmt.colheaders, 'LTpL_L1_r'));
LTpL_L2_r_ind = find(contains(lmt.colheaders, 'LTpL_L2_r'));
LTpL_L3_r_ind = find(contains(lmt.colheaders, 'LTpL_L3_r'));
LTpL_L4_r_ind = find(contains(lmt.colheaders, 'LTpL_L4_r'));
LTpL_L5_r_ind = find(contains(lmt.colheaders, 'LTpL_L5_r'));
IL_R5_r_ind = find(contains(lmt.colheaders, 'IL_R5_r'));
IL_R6_r_ind = find(contains(lmt.colheaders, 'IL_R6_r'));
IL_R7_r_ind = find(contains(lmt.colheaders, 'IL_R7_r'));
IL_R8_r_ind = find(contains(lmt.colheaders, 'IL_R8_r'));
IL_R9_r_ind = find(contains(lmt.colheaders, 'IL_R9_r'));
IL_R10_r_ind = find(contains(lmt.colheaders, 'IL_R10_r'));
IL_R11_r_ind = find(contains(lmt.colheaders, 'IL_R11_r'));
IL_R12_r_ind = find(contains(lmt.colheaders, 'IL_R12_r'));

LTpL_L1_l_ind = find(contains(lmt.colheaders, 'LTpL_L1_l'));
LTpL_L2_l_ind = find(contains(lmt.colheaders, 'LTpL_L2_l'));
LTpL_L3_l_ind = find(contains(lmt.colheaders, 'LTpL_L3_l'));
LTpL_L4_l_ind = find(contains(lmt.colheaders, 'LTpL_L4_l'));
LTpL_L5_l_ind = find(contains(lmt.colheaders, 'LTpL_L5_l'));
IL_R5_l_ind = find(contains(lmt.colheaders, 'IL_R5_l'));
IL_R6_l_ind = find(contains(lmt.colheaders, 'IL_R6_l'));
IL_R7_l_ind = find(contains(lmt.colheaders, 'IL_R7_l'));
IL_R8_l_ind = find(contains(lmt.colheaders, 'IL_R8_l'));
IL_R9_l_ind = find(contains(lmt.colheaders, 'IL_R9_l'));
IL_R10_l_ind = find(contains(lmt.colheaders, 'IL_R10_l'));
IL_R11_l_ind = find(contains(lmt.colheaders, 'IL_R11_l'));
IL_R12_l_ind = find(contains(lmt.colheaders, 'IL_R12_l'));

muscle_inds = [rect_abd_r_ind, rect_abd_l_ind, IL_L1_r_ind, IL_L1_l_ind, LTpT_T7_r_ind, LTpT_T7_l_ind, LTpL_L1_r_ind, LTpL_L1_l_ind];

IL_inds = [IL_L1_r_ind, IL_L2_r_ind, IL_L3_r_ind, IL_L4_r_ind, IL_L1_l_ind, IL_L2_l_ind, IL_L3_l_ind, IL_L4_l_ind];

IL_r_inds = [IL_L1_r_ind, IL_L2_r_ind, IL_L3_r_ind, IL_L4_r_ind];

IL_l_inds = [IL_L1_l_ind, IL_L2_l_ind, IL_L3_l_ind, IL_L4_l_ind];

LTpT_r_inds = [LTpT_T7_r_ind, LTpT_T8_r_ind, LTpT_T9_r_ind, LTpT_T10_r_ind, LTpT_T11_r_ind, LTpT_T12_r_ind, LTpT_R7_r_ind, LTpT_R8_r_ind, LTpT_R9_r_ind, LTpT_R10_r_ind, LTpT_R11_r_ind, LTpT_R12_r_ind];

LTpT_l_inds = [LTpT_T7_l_ind, LTpT_T8_l_ind, LTpT_T9_l_ind, LTpT_T10_l_ind, LTpT_T11_l_ind, LTpT_T12_l_ind, LTpT_R7_l_ind, LTpT_R8_l_ind, LTpT_R9_l_ind, LTpT_R10_l_ind, LTpT_R11_l_ind, LTpT_R12_l_ind];

LTpT_inds = [LTpT_T7_r_ind, LTpT_T8_r_ind, LTpT_T9_r_ind, LTpT_T10_r_ind, LTpT_T11_r_ind, LTpT_T12_r_ind, LTpT_R7_r_ind, LTpT_R8_r_ind, LTpT_R9_r_ind, LTpT_R10_r_ind, LTpT_R11_r_ind, LTpT_R12_r_ind, LTpT_T7_l_ind, LTpT_T8_l_ind, LTpT_T9_l_ind, LTpT_T10_l_ind, LTpT_T11_l_ind, LTpT_T12_l_ind, LTpT_R7_l_ind, LTpT_R8_l_ind, LTpT_R9_l_ind, LTpT_R10_l_ind, LTpT_R11_l_ind, LTpT_R12_l_ind];

LTpL_r_inds = [IL_R5_r_ind, IL_R6_r_ind, IL_R7_r_ind, IL_R8_r_ind, IL_R9_r_ind, IL_R10_r_ind, IL_R11_r_ind, IL_R12_r_ind, LTpL_L5_r_ind, LTpL_L4_r_ind, LTpL_L3_r_ind, LTpL_L2_r_ind, LTpL_L1_r_ind];

LTpL_l_inds = [IL_R5_l_ind, IL_R6_l_ind, IL_R7_l_ind, IL_R8_l_ind, IL_R9_l_ind, IL_R10_l_ind, IL_R11_l_ind, IL_R12_l_ind, LTpL_L5_l_ind, LTpL_L4_l_ind, LTpL_L3_l_ind, LTpL_L2_l_ind, LTpL_L1_l_ind];

LTpL_inds = [IL_R5_r_ind, IL_R6_r_ind, IL_R7_r_ind, IL_R8_r_ind, IL_R9_r_ind, IL_R10_r_ind, IL_R11_r_ind, IL_R12_r_ind, IL_R5_l_ind, IL_R6_l_ind, IL_R7_l_ind, IL_R8_l_ind, IL_R9_l_ind, IL_R10_l_ind, IL_R11_l_ind, IL_R12_l_ind, LTpL_L5_r_ind, LTpL_L4_r_ind, LTpL_L3_r_ind, LTpL_L2_r_ind, LTpL_L1_r_ind, LTpL_L5_l_ind, LTpL_L4_l_ind, LTpL_L3_l_ind, LTpL_L2_l_ind, LTpL_L1_l_ind];

%% find the muscle IDs of MA
rect_abd_r_ind = find(contains(ma.colheaders, 'rect_abd_r'));
rect_abd_l_ind = find(contains(ma.colheaders, 'rect_abd_l'));

IL_L1_r_ind = find(contains(ma.colheaders, 'IL_L1_r'));
IL_L2_r_ind = find(contains(ma.colheaders, 'IL_L2_r'));
IL_L3_r_ind = find(contains(ma.colheaders, 'IL_L3_r'));
IL_L4_r_ind = find(contains(ma.colheaders, 'IL_L4_r'));

IL_L1_l_ind = find(contains(ma.colheaders, 'IL_L1_l'));
IL_L2_l_ind = find(contains(ma.colheaders, 'IL_L2_l'));
IL_L3_l_ind = find(contains(ma.colheaders, 'IL_L3_l'));
IL_L4_l_ind = find(contains(ma.colheaders, 'IL_L4_l'));

LTpT_T7_r_ind = find(contains(ma.colheaders, 'LTpT_T7_r'));
LTpT_T7_l_ind = find(contains(ma.colheaders, 'LTpT_T7_l'));
LTpL_L1_r_ind = find(contains(ma.colheaders, 'LTpL_L1_r'));
LTpL_L1_l_ind = find(contains(ma.colheaders, 'LTpL_L1_l'));

ma_muscle_inds = [rect_abd_r_ind, rect_abd_l_ind, IL_L1_r_ind, IL_L1_l_ind, LTpT_T7_r_ind, LTpT_T7_l_ind, LTpL_L1_r_ind, LTpL_L1_l_ind];

%% Plotting muscle-tendon lengths

IL_r_LMT_sum = 0;
IL_l_LMT_sum = 0;
LTpT_r_LMT_sum = 0;
LTpT_l_LMT_sum = 0;
LTpL_r_LMT_sum = 0;
LTpL_l_LMT_sum = 0;

figure;
for i = 1:length(muscle_inds)
    muscle_ind = muscle_inds(i);
    subplot(4, 2, i);
    if i == 3
        for j = 1:length(IL_r_inds)
            IL_r_ind = IL_r_inds(j);
            plot(lmt.data(:, 1) - lmt.data(1, 1), lmt.data(:, IL_r_ind));
            hold on;
            IL_r_LMT_sum = IL_r_LMT_sum + lmt.data(:, IL_r_ind);
            if j == length(IL_r_inds)
                plot(lmt.data(:, 1) - lmt.data(1, 1), IL_r_LMT_sum / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 4
        for j = 1:length(IL_l_inds)
            IL_l_ind = IL_l_inds(j);
            plot(lmt.data(:, 1) - lmt.data(1, 1), lmt.data(:, IL_l_ind));
            hold on;
            IL_l_LMT_sum = IL_l_LMT_sum + lmt.data(:, IL_l_ind);
            if j == length(IL_l_inds)
                plot(lmt.data(:, 1) - lmt.data(1, 1), IL_l_LMT_sum / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 5
        for j = 1:length(LTpT_r_inds)
            LTpT_r_ind = LTpT_r_inds(j);
            plot(lmt.data(:, 1) - lmt.data(1, 1), lmt.data(:, LTpT_r_ind));
            hold on;
            LTpT_r_LMT_sum = LTpT_r_LMT_sum + lmt.data(:, LTpT_r_ind);
            if j == length(LTpT_r_inds)
                plot(lmt.data(:, 1) - lmt.data(1, 1), LTpT_r_LMT_sum / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 6
        for j = 1:length(LTpT_l_inds)
            LTpT_l_ind = LTpT_l_inds(j);
            plot(lmt.data(:, 1) - lmt.data(1, 1), lmt.data(:, LTpT_l_ind));
            hold on;
            LTpT_l_LMT_sum = LTpT_l_LMT_sum + lmt.data(:, LTpT_l_ind);
            if j == length(LTpT_l_inds)
                plot(lmt.data(:, 1) - lmt.data(1, 1), LTpT_l_LMT_sum / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 7
        for j = 1:length(LTpL_r_inds)
            LTpL_r_ind = LTpL_r_inds(j);
            plot(lmt.data(:, 1) - lmt.data(1, 1), lmt.data(:, LTpL_r_ind));
            hold on;
            LTpL_r_LMT_sum = LTpL_r_LMT_sum + lmt.data(:, LTpL_r_ind);
            if j == length(LTpL_r_inds)
                plot(lmt.data(:, 1) - lmt.data(1, 1), LTpL_r_LMT_sum / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 8
        for j = 1:length(LTpL_r_inds)
            LTpL_r_ind = LTpL_r_inds(j);
            plot(lmt.data(:, 1) - lmt.data(1, 1), lmt.data(:, LTpL_r_ind));
            hold on;
            LTpL_l_LMT_sum = LTpL_l_LMT_sum + lmt.data(:, LTpL_r_ind);
            if j == length(LTpL_r_inds)
                plot(lmt.data(:, 1) - lmt.data(1, 1), LTpL_l_LMT_sum / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    end

    title([lmt.colheaders{muscle_ind}]);
    ylabel('MuscleTendonLength');
end

total_muscle_num = length(IL_r_inds) + length(IL_l_inds) + length(LTpT_r_inds) + length(LTpT_l_inds) + length(LTpL_r_inds) + length(LTpL_l_inds);

Target_LMT = (IL_r_LMT_sum + IL_l_LMT_sum + LTpT_r_LMT_sum + LTpT_l_LMT_sum + LTpL_r_LMT_sum + LTpL_l_LMT_sum) / total_muscle_num;

%% Plotting muscle-tendon Forces

figure;
IL_r_force_SUM = 0;
IL_l_force_SUM = 0;
LTpT_r_force_SUM = 0;
LTpT_l_force_SUM = 0;
LTpL_r_force_SUM = 0;
LTpL_l_force_SUM = 0;
for i = 1:length(muscle_inds)
    muscle_ind = muscle_inds(i);
    subplot(4, 2, i);

    if i == 3
        for j = 1:length(IL_r_inds)
            IL_r_ind = IL_r_inds(j);
            plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, IL_r_ind));
            IL_r_force_SUM = IL_r_force_SUM + fmt.data(:, IL_r_ind);
            hold on;
            if j == length(IL_r_inds)
                plot(fmt.data(:, 1) - fmt.data(1, 1), IL_r_force_SUM, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 4
        for j = 1:length(IL_l_inds)
            IL_l_ind = IL_l_inds(j);
            plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, IL_l_ind));
            IL_l_force_SUM = IL_l_force_SUM + fmt.data(:, IL_l_ind);
            hold on;
            if j == length(IL_l_inds)
                plot(fmt.data(:, 1) - fmt.data(1, 1), IL_l_force_SUM, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 5
        for j = 1:length(LTpT_r_inds)
            LTpT_r_ind = LTpT_r_inds(j);
            plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, LTpT_r_ind));
            LTpT_r_force_SUM = LTpT_r_force_SUM + fmt.data(:, LTpT_r_ind);
            hold on;
            if j == length(LTpT_r_inds)
                plot(fmt.data(:, 1) - fmt.data(1, 1), LTpT_r_force_SUM, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 6
        for j = 1:length(LTpT_l_inds)
            LTpT_l_ind = LTpT_l_inds(j);
            plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, LTpT_l_ind));
            LTpT_l_force_SUM = LTpT_l_force_SUM + fmt.data(:, LTpT_l_ind);
            hold on;
            if j == length(LTpT_l_inds)
                plot(fmt.data(:, 1) - fmt.data(1, 1), LTpT_l_force_SUM, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 7
        for j = 1:length(LTpL_r_inds)
            LTpL_r_ind = LTpL_r_inds(j);
            plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, LTpL_r_ind));
            LTpL_r_force_SUM = LTpL_r_force_SUM + fmt.data(:, LTpL_r_ind);
            hold on;
            if j == length(LTpL_r_inds)
                plot(fmt.data(:, 1) - fmt.data(1, 1), LTpL_r_force_SUM, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 8
        for j = 1:length(LTpL_r_inds)
            LTpL_r_ind = LTpL_r_inds(j);
            plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, LTpL_r_ind));
            LTpL_l_force_SUM = LTpL_l_force_SUM + fmt.data(:, LTpL_r_ind);
            hold on;
            if j == length(LTpL_r_inds)
                plot(fmt.data(:, 1) - fmt.data(1, 1), LTpL_l_force_SUM, 'LineWidth', 2, 'Color', 'r');
            end
        end
    end

    title([fmt.colheaders{muscle_ind}]);
    ylabel('MuscleTendonForce');
end

%% Plotting fibre lengths

figure;
for i = 1:length(muscle_inds)
    muscle_ind = muscle_inds(i);
    subplot(4, 2, i);
    plot(lm.data(:, 1) - lm.data(1, 1), lm.data(:, muscle_ind));

    title([lm.colheaders{muscle_ind}]);
    ylabel('FibreLength');
end

%% Plotting tendon lengths

figure;
for i = 1:length(muscle_inds)
    muscle_ind = muscle_inds(i);
    subplot(4, 2, i);
    plot(lt.data(:, 1) - lt.data(1, 1), lt.data(:, muscle_ind));

    title([lt.colheaders{muscle_ind}]);
    ylabel('TendonLength');
end

%% Plotting fibre velocities

figure;
for i = 1:length(muscle_inds)
    muscle_ind = muscle_inds(i);
    subplot(4, 2, i);
    plot(vm.data(:, 1) - vm.data(1, 1), vm.data(:, muscle_ind));

    title([vm.colheaders{muscle_ind}]);
    ylabel('FibreVelocity');
end

figure;
plot(vm.data(:, 1) - vm.data(1, 1), vm.data(:, muscle_inds(1)));

%% Plotting moment arms around L5S1

figure;
IL_r_ma_SUM = 0;
IL_l_ma_SUM = 0;
LTpT_r_ma_SUM = 0;
LTpT_l_ma_SUM = 0;
LTpL_r_ma_SUM = 0;
LTpL_l_ma_SUM = 0;
for i = 1:length(ma_muscle_inds)
    muscle_ind = ma_muscle_inds(i);
    subplot(4, 2, i);

    if i == 3
        for j = 1:length(IL_r_inds)
            IL_r_ind = IL_r_inds(j);
            plot(ma.data(:, 1) - ma.data(1, 1), ma.data(:, IL_r_ind));
            IL_r_ma_SUM = IL_r_ma_SUM + ma.data(:, IL_r_ind);
            hold on;
            if j == length(IL_r_inds)
                plot(ma.data(:, 1) - ma.data(1, 1), IL_r_ma_SUM / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 4
        for j = 1:length(IL_l_inds)
            IL_l_ind = IL_l_inds(j);
            plot(ma.data(:, 1) - ma.data(1, 1), ma.data(:, IL_l_ind));
            IL_l_ma_SUM = IL_l_ma_SUM + ma.data(:, IL_l_ind);
            hold on;
            if j == length(IL_l_inds)
                plot(ma.data(:, 1) - ma.data(1, 1), IL_l_ma_SUM / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 5
        for j = 1:length(LTpT_r_inds)
            LTpT_r_ind = LTpT_r_inds(j);
            plot(ma.data(:, 1) - ma.data(1, 1), ma.data(:, LTpT_r_ind));
            LTpT_r_ma_SUM = LTpT_r_ma_SUM + ma.data(:, LTpT_r_ind);
            hold on;
            if j == length(LTpT_r_inds)
                plot(ma.data(:, 1) - ma.data(1, 1), LTpT_r_ma_SUM / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 6
        for j = 1:length(LTpT_l_inds)
            LTpT_l_ind = LTpT_l_inds(j);
            plot(ma.data(:, 1) - ma.data(1, 1), ma.data(:, LTpT_l_ind));
            LTpT_l_ma_SUM = LTpT_l_ma_SUM + ma.data(:, LTpT_l_ind);
            hold on;
            if j == length(LTpT_l_inds)
                plot(ma.data(:, 1) - ma.data(1, 1), LTpT_l_ma_SUM / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 7
        for j = 1:length(LTpL_r_inds)
            LTpL_r_ind = LTpL_r_inds(j);
            plot(ma.data(:, 1) - ma.data(1, 1), ma.data(:, LTpL_r_ind));
            LTpL_r_ma_SUM = LTpL_r_ma_SUM + ma.data(:, LTpL_r_ind);
            hold on;
            if j == length(LTpL_r_inds)
                plot(ma.data(:, 1) - ma.data(1, 1), LTpL_r_ma_SUM / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    elseif i == 8
        for j = 1:length(LTpL_r_inds)
            LTpL_r_ind = LTpL_r_inds(j);
            plot(ma.data(:, 1) - ma.data(1, 1), ma.data(:, LTpL_r_ind));
            LTpL_l_ma_SUM = LTpL_l_ma_SUM + ma.data(:, LTpL_r_ind);
            hold on;
            if j == length(LTpL_r_inds)
                plot(ma.data(:, 1) - ma.data(1, 1), LTpL_l_ma_SUM / j, 'LineWidth', 2, 'Color', 'r');
            end
        end
    end

    title([ma.colheaders{muscle_ind}]);
    ylabel('MomentArm');
end

%% torque
figure;

Torque_real = Tq.data(:, 2);
plot(Tq.data(:, 1) - Tq.data(1, 1), Tq.data(:, 2));

% title([ma.colheaders{'torque'}])
ylabel('L5S1_torque');

hold on;

IL_r_torque_SUM = 0;
IL_l_torque_SUM = 0;
LTpT_r_torque_SUM = 0;
LTpT_l_torque_SUM = 0;
LTpL_r_torque_SUM = 0;
LTpL_l_torque_SUM = 0;
torque_SUM = 0;

for i = 1:length(IL_r_inds)
    IL_r_ind = IL_r_inds(i);
    IL_r_torque_SUM = IL_r_torque_SUM + fmt.data(:, IL_r_ind) .* ma.data(:, IL_r_ind);
    plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, IL_r_ind) .* ma.data(:, IL_r_ind));
    hold on;
end
for i = 1:length(IL_l_inds)
    IL_l_ind = IL_l_inds(i);
    IL_l_torque_SUM = IL_l_torque_SUM + fmt.data(:, IL_l_ind) .* ma.data(:, IL_l_ind);
    plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, IL_l_ind) .* ma.data(:, IL_l_ind));
    hold on;
end
for i = 1:length(LTpT_r_inds)
    LTpT_r_ind = LTpT_r_inds(i);
    LTpT_r_torque_SUM = LTpT_r_torque_SUM + fmt.data(:, LTpT_r_ind) .* ma.data(:, LTpT_r_ind);
    plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, LTpT_r_ind) .* ma.data(:, LTpT_r_ind));
    hold on;
end
for i = 1:length(LTpT_l_inds)
    LTpT_l_ind = LTpT_l_inds(i);
    LTpT_l_torque_SUM = LTpT_l_torque_SUM + fmt.data(:, LTpT_l_ind) .* ma.data(:, LTpT_l_ind);
    plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, LTpT_l_ind) .* ma.data(:, LTpT_l_ind));
    hold on;
end
for i = 1:length(LTpL_r_inds)
    LTpL_r_ind = LTpL_r_inds(i);
    LTpL_r_torque_SUM = LTpL_r_torque_SUM + fmt.data(:, LTpL_r_ind) .* ma.data(:, LTpL_r_ind);
    plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, LTpL_r_ind) .* ma.data(:, LTpL_r_ind));
    hold on;
end
for i = 1:length(LTpL_l_inds)
    LTpL_l_ind = LTpL_l_inds(i);
    LTpL_l_torque_SUM = LTpL_l_torque_SUM + fmt.data(:, LTpL_l_ind) .* ma.data(:, LTpL_l_ind);
    plot(fmt.data(:, 1) - fmt.data(1, 1), fmt.data(:, LTpL_l_ind) .* ma.data(:, LTpL_l_ind));
    hold on;
end

torque_SUM = IL_r_torque_SUM + LTpT_r_torque_SUM + LTpL_r_torque_SUM + IL_l_torque_SUM + LTpT_l_torque_SUM + LTpL_l_torque_SUM;

plot(fmt.data(:, 1) - fmt.data(1, 1), torque_SUM);

activation_W(:, 1) = IL_r_torque_SUM ./ torque_SUM;
activation_W(:, 2) = IL_l_torque_SUM ./ torque_SUM;
activation_W(:, 3) = LTpT_r_torque_SUM ./ torque_SUM;
activation_W(:, 4) = LTpT_l_torque_SUM ./ torque_SUM;
activation_W(:, 5) = LTpL_r_torque_SUM ./ torque_SUM;
activation_W(:, 6) = LTpL_l_torque_SUM ./ torque_SUM;

%% the data sychronize function flag 1, flag 2, flag 3, and flag 4

for z = 1:length(ikData.data(:, 1))
    if lmt.data(1, 1) == ikData.data(z, 1)
        flag1 = z;
    end
    if lmt.data(end, 1) == ikData.data(z, 1)
        flag2 = z;
    end
end

for w = 1:length(Ac.data(:, 1))
    if lmt.data(5, 1) - Ac.data(w, 1) >= 0
        flag3 = w;
    end
    if lmt.data(end, 1) - Ac.data(w, 1) >= 0
        flag4 = w;
    end
end

%% Compute the Muscle tendon force

ikdata_all_comp = interp1(ikData.data(flag1:flag2, 1), ikData.data(flag1:flag2, 2), lmt.data(:, 1));

offset = -10 * ikdata_all_comp(1);

theta = offset + pi / 2 + 10 * ikdata_all_comp;

theta(1) = theta(3);
theta(2) = theta(3);

% L2 = 0.17280;
% L3 = 0.006;
% R = 0.052;
%
% [LMTT, MA] = simplify_pulley_model(L2, L3, R, theta*180/pi);

LMTT(2) = 0.164700000000000;
% 0.164700000000000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EMGGG = interp1(EMGfilt.data(4705:22881,1), EMGfilt.data(4705:22881,4), lmt.data(206:758,1));
% activation = interp1(Ac.data(1:840,1), Ac.data(1:840,muscle_inds(3)), lmt.data(1:758,1));
SUM_activation = 0;

activation = 0;

% for i = 4:9
%     activations(:,i-3) = Ac.data(:,muscle_ind);
%     %activations(:,i-3) = interp1(EMGfilt.data(flag3:flag4,1), EMGfilt.data(flag3:flag4,i), lmt.data(:,1));
%     SUM_activation = SUM_activation + activations(:,i-3);
% end

for i = 3:length(muscle_inds)
    muscle_ind = muscle_inds(i);
    activations(:, i - 2) = interp1(Ac.data(flag3:flag4, 1), Ac.data(flag3:flag4, muscle_ind), lmt.data(:, 1));
    SUM_activation = SUM_activation + activations(:, i - 2);
end

for i = 1:6
    activation = activation + activations(:, i) .* activation_W(:, i);
end

% NaN filter
for j = 1:length(activation) %%%%
    if isnan(activation(j))
        activation(j) = 0.015;
    end
end

t_start = lmt.data(153, 1); % 起始时间%143
t_end = lmt.data(267, 1); % 结束时间248

t_interp = (t_start:0.002:t_end)';

theta = interp1(lmt.data(:, 1), theta, t_interp);

Torque_real = interp1(lmt.data(:, 1), Torque_real, t_interp);

activation = interp1(lmt.data(:, 1), activation, t_interp);

Smooth_activation = smoothdata(activation, "gaussian", 5);

MTF_real = fmt.data(:, IL_L1_r_ind);

LoM_nominal = 0.056058545564031212;

LsT = 0.12706603661180407;

damping = 0.5;

alpha0 = 0.24086; % obtain from the paper

FoM = 962.3;
% 1924.6

Ts = 0.002;

L2 = 0.17;

L3 = 0.005;

R = 0.05;

% model_torque = Torque_function(theta, activation, damping, alpha0, Ts, L2, L3, R, LoM_nominal, LsT, FoM);
%
% %% using the least square method to obtain the optimal parameters
%
% objective = @(coeffs) norm(Torque_function(theta, activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6)) - Torque_real);
%
% %objective_value = Torque_fitness_function(Torque_real(206:758), theta, activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6));
%
%
%
% initial_guess = [0.1728 0.005 0.052 0.056058545564031212 0.12706603661180407 4000];
% %initial_guess = [0.1728 0.005 0.052 0.056058545564031212 0.12706603661180407 7000]0.11 0.0051 0.033 0.056058545564031212 0.12706603661180407 6000;
% lb = [0.1 0.002 0.02 0 0 3000];
% ub = [1 0.1 0.2 0.2 1 15000];
% A = [-1 0 0 1 1 0];
% b = 0.03;
% Aeq = [];
% beq = [];
% nonlcon = [];
% %options = optimoptions('lsqnonlin', 'Display', 'iter','Algorithm','trust-region-reflective','TolX', 1e-14);
%
% options = optimoptions('ga','MaxGenerations', 1500 ,'PlotFcn', @gaplotbestf);

% coeffs = ga(objective, 6, A, b, Aeq, beq, lb, ub,nonlcon , options);
% coeffs  = lsqnonlin(objective,initial_guess,lb,ub,options);
% coeffs  = lsqnonlin(objective, initial_guess, lb, ub, A, b, Aeq, beq, @LMT_Constraint, options);
%%

coeffs = [0.251891307986178 0.00819647269406587 0.0996898296899076  0.0364835264780155  0.207132716708542   3794.97267881575    3.82557174688585    7.12655438615596    -78.2474304632419   -81.1892832678883];

[LMT, BE] = simplify_pulley_model(coeffs(1), coeffs(2), coeffs(3), theta);

Model_torque = Torque_function(Smooth_activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10));

% Model_MTF = MTU_function(LMT, theta, Smooth_activation, damping, alpha0, Ts, coeffs(1), coeffs(2), coeffs(3), coeffs(4), coeffs(5), coeffs(6), coeffs(7), coeffs(8), coeffs(9), coeffs(10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the joint torque obtain by one rebuild muscle

figure;

plot(t_interp - t_start, Torque_real);

hold on;

plot(t_interp - t_start, Model_torque);

ylabel('L5S1_torque');

legend('Real-torque', 'Model-torque');
title('Real-torque, Model-torque');

figure;
plot(t_interp - t_start, activation);
hold on;
% plot(fmt.data(206:758,1) - fmt.data(206,1), EMGGG);
% legend('Average activation', 'EMGfilt-IL-r');
title('Average activation, EMGfilt-IL-r');

%% plot the model_LMT and target_LMT
figure;

plot(t_interp - t_start, interp1(lmt.data(:, 1), Target_LMT, t_interp));

hold on;

plot(t_interp - t_start, LMT);

ylabel('Muscle tendon length');

legend('Target/mean-LMT', 'Model-LMT');
title('Target/mean-LMT, Model-LMT');

%% plot the moment of arm

figure;
plot(t_interp - t_start, BE);
title('MA');

%% plot the all the outputs together
% torque
R = corrcoef(Model_torque, Torque_real);
correlation(1) = R(1, 2);
figure('DefaultAxesFontSize', 24, 'DefaultTextFontSize', 24);
subplot(2, 1, 1);
plot(t_interp - t_start, Torque_real, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3.5);

hold on;

plot(t_interp - t_start, Model_torque, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5);
xlabel('Time (s)');
ylabel('Torque(N.M)');

Model_torque_max = max(Model_torque); % rescale the profile
ylim([0, Model_torque_max] * 1.3);
legend('CEINMS-Torque', 'Model-Torque', 'FontSize', 20);
title('(b)                   ');
legend boxoff;
% % MTF
% subplot(2,1,2);
% plot(t_interp - t_start, Model_MTF, 'LineWidth', 2);
% xlabel('Time (s)');
% ylabel('Force(N)')
%
% Model_MTF_max = max(Model_MTF); % rescale the profile
% ylim([0, Model_MTF_max] * 1.1);
%
% legend('MTF','FontSize', 15);
% title('(c) musculotendon force generate by the model');

% LMT

R = corrcoef(LMT, interp1(lmt.data(:, 1), Target_LMT, t_interp));
correlation(2) = R(1, 2);
subplot(2, 1, 2);
plot(t_interp - t_start, interp1(lmt.data(:, 1), Target_LMT, t_interp), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3.5);

hold on;

plot(t_interp - t_start, LMT, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 3.5);

xlabel('Time (s)');
ylabel('L^{MT}(m)');
LMT_max = max(LMT) - min(Target_LMT); % rescale the profile
LMT_min = min(Target_LMT); % rescale the profile
ylim([LMT_min, LMT_min + LMT_max * 1.3]);
legend('CEINMS-L^{MT}', 'Model-L^{MT}', 'FontSize', 20);
title('(c)                      ');
legend boxoff;
% % moment arm
%
% subplot(2,2,4);
% plot(t_interp - t_start, BE, 'LineWidth', 2);
% xlabel('Time (s)');
%
% BE_max = max(BE) - min(BE); % rescale the profile
% BE_min = min(BE); % rescale the profile
% ylim([BE_min, BE_min + BE_max * 1.1]);
%
% legend('MA','FontSize', 15);
% ylabel('MA(m)')
% title('(d) Moment arm of the model');

% 调整子图的位置和宽度，使其保持一致
subplot_handles = findall(gcf, 'Type', 'Axes');
for i = 1:length(subplot_handles)
    subplot_handles(i).Position(3) = 0.3; % 设置子图的宽度
end
