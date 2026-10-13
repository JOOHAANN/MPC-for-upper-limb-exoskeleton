function [Reference, Scale] = Dynamic_kine_reference(childFolder, startpoint, endpoint)
    % DYNAMIC_KINE_REFERENCE Extract synchronized reference signals from experimental recordings.
    % Use prepare_workflow to select matching code and input data.

    parentFolder = exo_data_root();

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
        if i == 3
            for j = 1:length(IL_r_inds)
                IL_r_ind = IL_r_inds(j);
                IL_r_LMT_sum = IL_r_LMT_sum + lmt.data(:, IL_r_ind);
            end
        elseif i == 4
            for j = 1:length(IL_l_inds)
                IL_l_ind = IL_l_inds(j);
                IL_l_LMT_sum = IL_l_LMT_sum + lmt.data(:, IL_l_ind);
            end
        elseif i == 5
            for j = 1:length(LTpT_r_inds)
                LTpT_r_ind = LTpT_r_inds(j);
                LTpT_r_LMT_sum = LTpT_r_LMT_sum + lmt.data(:, LTpT_r_ind);
            end
        elseif i == 6
            for j = 1:length(LTpT_l_inds)
                LTpT_l_ind = LTpT_l_inds(j);
                LTpT_l_LMT_sum = LTpT_l_LMT_sum + lmt.data(:, LTpT_l_ind);
            end
        elseif i == 7
            for j = 1:length(LTpL_r_inds)
                LTpL_r_ind = LTpL_r_inds(j);
                LTpL_r_LMT_sum = LTpL_r_LMT_sum + lmt.data(:, LTpL_r_ind);
            end
        elseif i == 8
            for j = 1:length(LTpL_r_inds)
                LTpL_r_ind = LTpL_r_inds(j);
                LTpL_l_LMT_sum = LTpL_l_LMT_sum + lmt.data(:, LTpL_r_ind);
            end
        end
    end

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

        if i == 3
            for j = 1:length(IL_r_inds)
                IL_r_ind = IL_r_inds(j);
                IL_r_force_SUM = IL_r_force_SUM + fmt.data(:, IL_r_ind);
            end
        elseif i == 4
            for j = 1:length(IL_l_inds)
                IL_l_ind = IL_l_inds(j);
                IL_l_force_SUM = IL_l_force_SUM + fmt.data(:, IL_l_ind);
            end
        elseif i == 5
            for j = 1:length(LTpT_r_inds)
                LTpT_r_ind = LTpT_r_inds(j);
                LTpT_r_force_SUM = LTpT_r_force_SUM + fmt.data(:, LTpT_r_ind);
            end
        elseif i == 6
            for j = 1:length(LTpT_l_inds)
                LTpT_l_ind = LTpT_l_inds(j);
                LTpT_l_force_SUM = LTpT_l_force_SUM + fmt.data(:, LTpT_l_ind);
            end
        elseif i == 7
            for j = 1:length(LTpL_r_inds)
                LTpL_r_ind = LTpL_r_inds(j);
                LTpL_r_force_SUM = LTpL_r_force_SUM + fmt.data(:, LTpL_r_ind);
            end
        elseif i == 8
            for j = 1:length(LTpL_r_inds)
                LTpL_r_ind = LTpL_r_inds(j);
                LTpL_l_force_SUM = LTpL_l_force_SUM + fmt.data(:, LTpL_r_ind);
            end
        end
    end

    total_muscle_num = length(IL_r_inds) + length(IL_l_inds) + length(LTpT_r_inds) + length(LTpT_l_inds) + length(LTpL_r_inds) + length(LTpL_l_inds);

    Target_LMT = (IL_r_LMT_sum + IL_l_LMT_sum + LTpT_r_LMT_sum + LTpT_l_LMT_sum + LTpL_r_LMT_sum + LTpL_l_LMT_sum) / total_muscle_num;

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

        if i == 3
            for j = 1:length(IL_r_inds)
                IL_r_ind = IL_r_inds(j);
                IL_r_ma_SUM = IL_r_ma_SUM + ma.data(:, IL_r_ind);
            end
        elseif i == 4
            for j = 1:length(IL_l_inds)
                IL_l_ind = IL_l_inds(j);
                IL_l_ma_SUM = IL_l_ma_SUM + ma.data(:, IL_l_ind);
            end
        elseif i == 5
            for j = 1:length(LTpT_r_inds)
                LTpT_r_ind = LTpT_r_inds(j);
                LTpT_r_ma_SUM = LTpT_r_ma_SUM + ma.data(:, LTpT_r_ind);
            end
        elseif i == 6
            for j = 1:length(LTpT_l_inds)
                LTpT_l_ind = LTpT_l_inds(j);
                LTpT_l_ma_SUM = LTpT_l_ma_SUM + ma.data(:, LTpT_l_ind);
            end
        elseif i == 7
            for j = 1:length(LTpL_r_inds)
                LTpL_r_ind = LTpL_r_inds(j);
                LTpL_r_ma_SUM = LTpL_r_ma_SUM + ma.data(:, LTpL_r_ind);
            end
        elseif i == 8
            for j = 1:length(LTpL_r_inds)
                LTpL_r_ind = LTpL_r_inds(j);
                LTpL_l_ma_SUM = LTpL_l_ma_SUM + ma.data(:, LTpL_r_ind);
            end
        end
    end

    IL_r_torque_SUM = 0;
    IL_l_torque_SUM = 0;
    LTpT_r_torque_SUM = 0;
    LTpT_l_torque_SUM = 0;
    LTpL_r_torque_SUM = 0;
    LTpL_l_torque_SUM = 0;
    torque_SUM = 0;

    Torque_real = Tq.data(:, 2);

    Scale(:, 1) = Tq.data(:, 2);

    for i = 1:length(IL_r_inds)
        IL_r_ind = IL_r_inds(i);
        IL_r_torque_SUM = IL_r_torque_SUM + fmt.data(:, IL_r_ind) .* ma.data(:, IL_r_ind);
    end
    for i = 1:length(IL_l_inds)
        IL_l_ind = IL_l_inds(i);
        IL_l_torque_SUM = IL_l_torque_SUM + fmt.data(:, IL_l_ind) .* ma.data(:, IL_l_ind);
    end
    for i = 1:length(LTpT_r_inds)
        LTpT_r_ind = LTpT_r_inds(i);
        LTpT_r_torque_SUM = LTpT_r_torque_SUM + fmt.data(:, LTpT_r_ind) .* ma.data(:, LTpT_r_ind);
    end
    for i = 1:length(LTpT_l_inds)
        LTpT_l_ind = LTpT_l_inds(i);
        LTpT_l_torque_SUM = LTpT_l_torque_SUM + fmt.data(:, LTpT_l_ind) .* ma.data(:, LTpT_l_ind);
    end
    for i = 1:length(LTpL_r_inds)
        LTpL_r_ind = LTpL_r_inds(i);
        LTpL_r_torque_SUM = LTpL_r_torque_SUM + fmt.data(:, LTpL_r_ind) .* ma.data(:, LTpL_r_ind);
    end
    for i = 1:length(LTpL_l_inds)
        LTpL_l_ind = LTpL_l_inds(i);
        LTpL_l_torque_SUM = LTpL_l_torque_SUM + fmt.data(:, LTpL_l_ind) .* ma.data(:, LTpL_l_ind);
    end

    torque_SUM = IL_r_torque_SUM + LTpT_r_torque_SUM + LTpL_r_torque_SUM + IL_l_torque_SUM + LTpT_l_torque_SUM + LTpL_l_torque_SUM;

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

    %% normlize the theta from 30 to 84
    ikdata_all_comp = interp1(ikData.data(flag1:flag2, 1), ikData.data(flag1:flag2, 2), lmt.data(:, 1));

    offset = max(ikdata_all_comp) + 0.02; % 0.01

    theta = pi / 2 + 10 * (ikdata_all_comp - offset); % 7

    theta(1) = theta(3);
    theta(2) = theta(3);
    %% normlize the theta from 30 to 84

    Scale(:, 2) = theta; % theta scale to the lmt, fmt, .....

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

    %% changing the sampling freq

    Ts = 0.001;

    t_start = lmt.data(startpoint, 1); % 起始时间177&119&154&123
    t_end = lmt.data(endpoint, 1); % 结束时间219&269&264&295

    t_interp = (t_start:Ts:t_end)';

    theta_raw = theta;

    theta = interp1(lmt.data(:, 1), theta, t_interp);

    theta = smoothdata(theta, "gaussian", 1500);

    Torque_real = interp1(lmt.data(:, 1), Torque_real, t_interp);

    activation = interp1(lmt.data(:, 1), activation, t_interp);

    Target_LMT = interp1(lmt.data(:, 1), Target_LMT, t_interp);

    Smooth_activation = smoothdata(activation, "gaussian", 5);

    %% plot the joint torque obtain by one rebuild muscle

    figure;

    plot(t_interp, Torque_real);

    ylabel('L5S1_torque');

    legend('Real-torque');
    title('Real-torque');

    figure;
    plot(t_interp, activation);

    title('Average activation, EMGfilt-IL-r');

    figure;
    plot(theta_raw);
    title('Raw theta');

    %% plot the model_LMT and target_LMT
    figure;

    plot(t_interp, Target_LMT);
    ylabel('Muscle tendon length');

    legend('Target-LMT');
    title('Target-LMT');

    Reference(:, 1) = Torque_real;
    Reference(:, 2) = Smooth_activation;
    Reference(:, 3) = theta;
    Reference(:, 4) = Target_LMT;

end
