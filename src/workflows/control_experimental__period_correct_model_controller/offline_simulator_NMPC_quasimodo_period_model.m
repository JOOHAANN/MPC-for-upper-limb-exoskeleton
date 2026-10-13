% OFFLINE_SIMULATOR_NMPC_QUASIMODO_PERIOD_MODEL Run the offline NMPC simulation for the selected activation source.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;

%% this is the offline NMPC simulator using in Quasimodo

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

LMTT(2) = 0.164700000000000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUM_activation = 0;

activation = 0;

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

t_start = lmt.data(153, 1); % start time %143
t_end = lmt.data(267, 1); % end time 248

t_interp = (t_start:0.002:t_end)';

theta = interp1(lmt.data(:, 1), theta, t_interp);

Torque_real = interp1(lmt.data(:, 1), Torque_real, t_interp);

activation = interp1(lmt.data(:, 1), activation, t_interp);

Smooth_activation = smoothdata(activation, "gaussian", 5);

%%%%%%%%% start with designing the NMPC %%%%%%%%%%%%%%
%% scale muscle model

coeffs = [0.254236968598434 0.00477342247054632 0.129274795804972   0.184459724281561   0.0951457946485656  5246.89114150244    129.977390953028    1.55323492490645    -295.243685003012   373.521985787704];

L2 = coeffs(1);
L3 = coeffs(2);
R = coeffs(3);
LoM_nominal = coeffs(4);
LsT = coeffs(5);
FoM = coeffs(6);
I = coeffs(7);
M = coeffs(8);
a = coeffs(9);
b = coeffs(10);

%% real-time NMPC

%% Preallocations

NumberOfData = length(Smooth_activation);

theta = zeros(NumberOfData, 1)';
dtheta = zeros(NumberOfData, 1)';
ddtheta = zeros(NumberOfData, 1)';

LoM = zeros(NumberOfData, 1)';
alpha = zeros(NumberOfData, 1)';
dalpha = zeros(NumberOfData, 1)';
V_max = zeros(NumberOfData, 1)';
LMT = zeros(NumberOfData, 1)';
MA = zeros(NumberOfData, 1)';
LM = zeros(NumberOfData, 1)';
dLM = zeros(NumberOfData, 1)';
LT = zeros(NumberOfData, 1)';
dLMT = zeros(NumberOfData, 1)';
k_t = zeros(NumberOfData, 1)';
stiffness = zeros(NumberOfData, 1)';
F_L = zeros(NumberOfData, 1)';
F_PE = zeros(NumberOfData, 1)';

t = zeros(NumberOfData, 1)';

MTF = zeros(NumberOfData + 1, 1)';
u_Tac = zeros(NumberOfData + 1, 1)';
u_Tac_coeff = zeros(NumberOfData + 1, 1)';
dFt = zeros(NumberOfData, 1)';
% Tac = zeros(NumberOfData,1)';
strain = zeros(NumberOfData, 1)';

%% initial the state

PercentageChange = 0.15;

Ts = 0.002;

g = 9.81;

kt_spline = [6.07450471822638e+42   -1.17093401235721e+42   9.78861307802160e+40    -4.68161016668128e+39   1.41111335336212e+38    -2.75023728531764e+36   3.32535388381590e+34    -1.97635568280683e+32   -4.74852712911846e+29   1.68812705982292e+28    -9.30681300424510e+25   -1.69227796153188e+23   3.64048821058250e+21    -7.76789223923884e+18   -4.91004587706162e+16   201470178531589 239829087287.116    -1720265398.11745   -15739.3275509684   6708.90359529906    3.90612013684262];

damping = 0.5;

alpha0 = 0.24086;

w = LoM_nominal * sin(alpha0);

InitialConditions = [10, 1.6, 0];

MTF(1) = InitialConditions(1); % state x1

theta(1) = InitialConditions(2); % state x2

dtheta(1) = InitialConditions(3); % state x3

[LMT(1), MA(1)] = simplify_pulley_model(L2, L3, R, theta(1));

u_Tac(1) = 0;

for k = 1:(NumberOfData - 1)

    %% compute LM for time k
    activation(k) = activation(k) * (1 - u_Tac_coeff(k));

    if k == 1
        k_t(k) = 10 * FoM / LsT;
    else
        strain(k) = (LT(k - 1) - LsT) / LsT;
        if strain(k) > 0.0199
            strain(k) = 0.0199;
            stiffness(k) = polyval(kt_spline, strain(k));
        else
            stiffness(k) = polyval(kt_spline, strain(k));
        end
        k_t(k) = stiffness(k) * FoM / LsT;
    end

    LoM(k) = LoM_nominal * (1 + PercentageChange * (1 - activation(k)));

    V_max(k) = 10 * LoM(k);

    LT(k) = LsT + MTF(k) / k_t(k);                                                % tendon length
    alpha(k) = acos(sqrt((LMT(k) - LT(k))^2 / (w^2 + (LMT(k) - LT(k))^2)));         % compute the pennation angles for the time k
    LM(k) = (LMT(k) - LT(k)) / cos(alpha(k));                                         % muscle length

    %% compute dLM

    F_L(k) = exp(-17.33 * (abs((LM(k) / LoM(k)) - 1)^3)); % equ (42)
    F_PE(k) = FoM * 3 * 10000 * exp(6 * ((LM(k) / LoM(k)) - 3.2)); % equ (12)
    % F_V(k)=(MTF(k) / cos(alpha(k))-F_PE(k)) / (activation(k) * FoM * exp(-17.33 * abs((LM(k) / LoM(k)) - 1)^3));%solve directly

    %% 1. directive method inverse solving

    %% 2. regression method

    e = [2.2 1];                                                             %  Fitted exactly to the linear part of Fv (k2 in the "Testing_Fv" file
    Z1 = MTF(k) / (FoM * cos(alpha(k))) - F_PE(k) / FoM - activation(k) * F_L(k) * e(2);
    Z2 = activation(k) * F_L(k) * e(1) / (V_max(k)) + damping / V_max(k);
    % dLM(k)=Vw_max * LoM * (0.5+exp(-3*(F_V(k)-0.02))+0.1*exp(12*(F_V(k)-1.5)));

    dLM(k) = Z1 / Z2;

    dalpha(k) = -dLM(k) / LM(k) * tan(alpha(k));

    u_Tac = MTF(k) * MA(k) * u_Tac_coeff(k) / (1 - u_Tac_coeff(k));

    ddtheta(k) = (MTF(k) * MA(k) + u_Tac + a * theta(k)^2 + b * theta(k) - M * cos(theta(k))) / I;

    %% integration and compute Ft

    theta(k + 1) = theta(k) + dtheta(k) * Ts;
    dtheta(k + 1) = dtheta(k) + ddtheta(k) * Ts;
    % LMT(k+1) = sqrt(L3^2 + L2^2 + 2*L2*L3*cos(theta(k+1)));
    [LMT(k + 1), MA(k + 1)] = simplify_pulley_model(L2, L3, R, theta(k + 1));
    dLMT(k + 1) = L2 * L3 * sin(theta(k + 1)) * dtheta(k + 1);

    dFt(k) = k_t(k) * (dLMT(k) + LM(k) * sin(alpha(k)) * dalpha(k) - dLM(k) * cos(alpha(k))); % handnotes equ (37) and (16)

    MTF(k + 1) = MTF(k) + dFt(k) * Ts;

    % LMT(k+1) = LMT(k) + dLMT(k)*Ts;

    if MTF(k + 1) < 0
        MTF(k + 1) = 0;
    end

    if isnan(MTF(k + 1))
        MTF(k + 1) = MTF(k);
    end

    if MTF(k + 1) > 7000
        MTF(k + 1) = MTF(k);
    end

    t(k + 1) = k * Ts;

    torque = MTF(k + 1) * MA(k + 1);
    %% optimization by cost function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NMPC

    if torque > 80
        Q = 100; % state weight
        RR = 1; % control weight
        InitialConditions = [MTF(k + 1), theta(k + 1), dtheta(k + 1)]; % initial state
        r = 110; % ref

        N = 50;

        cost_function = @(z) mpc_cost_function(z, InitialConditions, N, Q, RR, r, damping, alpha0, Ts, activation(k:k + N), coeffs);
        z0 = 50 * ones(N, 1);
        lb = -10;
        ub = 300;
        % using fmincon
        options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'UseParallel', true);
        % using GA
        % options = optimoptions('ga', 'Display', 'off', 'UseParallel', true);

        % using fmincon
        z = fmincon(cost_function, z0, [], [], [], [], 0 * ones(N, 1), 0.8 * ones(N, 1), [], options); % @(z) state_constraints(z, InitialConditions, N, Ts, 0, activation, coeffs, k_t(k))
        % using GA
        % z = ga(cost_function, N, [], [], [], [], lb, ub, [], options);

        % extract the first control input
        z_initials = z;
        u_Tac_coeff(k + 1) = z(1);

        % % apply the control input
        % [MTF_next, theta_next, dtheta_next] = nonlinear_system(k_t(k), InitialConditions, activation(k), u_Tac, Ts, coeffs);
        %
        x = [MTF(k + 1) * MA(k + 1), theta(k + 1), dtheta(k + 1)];
        % saving the result
        X(:, k) = x;
        U(k) = u_Tac_coeff(k + 1);
        Y(k) = x(1); % output the state x1
    else
        u_Tac_coeff(k + 1) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
model_torque = MTF(1:end - 1)' .* MA';
