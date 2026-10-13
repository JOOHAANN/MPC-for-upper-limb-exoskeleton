% NNMF_COMPONENT Estimate muscle synergies and their speed/load dependence from EMG.
% Use prepare_workflow to select matching code and input data.
clc;
clear;
close all;

Ts = 0.002;
fs = 1 / Ts;

set(0, 'DefaultFigureWindowStyle', 'docked');

%%%%%%%%%%%% Step 1%%%%%%%%%% define the path of the data folder
parentFolder = exo_data_root();

% set the path as default
baseFolder = '240528Pilot';
subFolder1 = 'CalibrationData';
movementType = 'stoop';

weight = [5, 75, 10, 125, 15]; % weight
speed = [25, 35, 45, 55, 65]; % speed

% Number of files to import
numFiles = 25;
EMGfilts = struct('data', cell(1, numFiles));
EMGParts = struct('data', cell(1, numFiles));
lmt = struct('data', cell(1, numFiles));
childFolders = cell(1, 25);
envelopeEMG = struct('data', cell(1, numFiles));
EMG_matrix = struct('data', cell(1, numFiles));
W = struct('data', cell(1, numFiles));
W_norm = struct('data', cell(1, numFiles));
H = struct('data', cell(1, numFiles));
xq = struct('data', cell(1, numFiles));
INW = struct('data', cell(1, numFiles));
Weightings = struct('data', cell(1, numFiles));
Increment_Speed = struct('data', cell(1, numFiles));
Increment_Weight = struct('data', cell(1, numFiles));

%%%%%%%%%%%% Step 2%%%%%%%%%% define the segment of the EMG signals
S = [11797, 17280; 10364, 14324; 2808, 5021; 4620, 6669; 3486, 5154
     4952, 9416; 5727, 8919; 3453, 6377; 4075, 6350; 3523, 5447
     5420, 9558; 3864, 7572; 5330, 8588; 2368, 3831; 3723, 5388
     3515, 7448; 2416, 5825; 3825, 7001; 3945, 6189; 2224, 3532
     5774, 9452; 4358, 7197; 4552, 7392; 1878, 3470; 4468, 6102];
%%%%%%%%%%%% Step 2%%%%%%%%%% define the segment of the EMG signals

i = 1;

for is = 1:5
    for js = 1:5
        subFolder2 = sprintf('dynamics_%s_%dkg_%d', movementType, weight(is), speed(js));
        childFolders{i} = fullfile(baseFolder, subFolder1, subFolder2); % combine the variables as total path
        % Preallocate a struct array to store the imported data
        EMGfilts(i).data = importdata(exo_data_file(parentFolder, childFolders{i}, 'emgFilt.sto'));
        for ip = 1:6
            EMGParts(i).data(:, ip) = EMGfilts(i).data.data(S(i, 1):S(i, 2), ip + 3);
        end
        EMGParts(i).data(:, 7) = EMGfilts(i).data.data(S(i, 1):S(i, 2), 1);
        lmt(i).data = importdata(exo_data_file(parentFolder, childFolders{i}, 'lmt.sto'));
        i = i + 1;
    end
end

lowpassCutoff = 3;  % lowpass filter cutoff freq
[d, c] = butter(4, lowpassCutoff / (fs / 2), 'low');

%% the data sychronize function flag 1, flag 2, flag 3, and flag 4

for ig = 1:numFiles

    StartP(ig) = EMGParts(ig).data(1, 7);

    EndP(ig) = EMGParts(ig).data(end, 7);

    xq(ig).data = linspace(StartP(ig), EndP(ig), 200);

    for iE = 1:6
        envelopeEMG(ig).data(:, iE) = filtfilt(d, c, EMGParts(ig).data(:, iE));
        % envelopeEMG(ig).data(:,iE) = smoothdata(EMGParts(ig).data(:,iE), "gaussian", 5);
        EMG_matrix(ig).data(:, iE) = interp1(EMGParts(ig).data(:, 7), envelopeEMG(ig).data(:, iE), xq(ig).data(1:200));
        % NaN filter
        for j = 1:200 %%%%
            if isnan(EMG_matrix(ig).data(j, iE))
                EMG_matrix(ig).data(j, iE) = EMG_matrix(ig).data(j - 1, iE);
            end
        end
    end

end

filename_EMG = 'EMG_matrix.mat';
save(filename_EMG, 'EMG_matrix');

k = 1:200;

SUM_W = 0;

figure;
for ix = 1:numFiles
    [W(ix).data, H(ix).data] = nnmf(EMG_matrix(ix).data, 4);
    % for ic = 1:4
    % W_norm(ix).data(:,ic) = normalize(W(ix).data(:,ic),'range');
    % end
    [W_norm_data, settings] = mapminmax(W(ix).data', 0, 1);
    W_norm(ix).data = W_norm_data';
    % W_norm(ix).data = rescale(W(ix).data);

    SUM_W =  SUM_W + W_norm(ix).data;
    if ix == 25
        Baselines = SUM_W ./ 25;
    end
    INW(ix).data = pinv(W_norm(ix).data);
    Weightings(ix).data = INW(ix).data * EMG_matrix(ix).data;
end

%% save the Weighting_baseline and XP_baseline data to the a file
baseline_data = Weightings(13).data;

filename_W_baseline = 'baseline_data.mat';
filename_XP_baseline = 'Baselines.mat';

% save the data in the baseline_file
save(filename_W_baseline, 'baseline_data');
save(filename_XP_baseline, 'Baselines');

%% the speed and weight-dependent increments

i = 1;

for ic = 1:5
    for jc = 1:5
        if jc == 3
            for il = 1:5
                Increment_Speed(i + il - 3).data = Weightings(i + il - 3).data - Weightings(i).data;
            end
        end
        if ic == 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% original value is ic = 3 ;(i+5*(il-1))
            for il = 1:5
                Increment_Weight(i + 5 * (il - 3)).data = Weightings(i + 5 * (il - 3)).data - Weightings(i).data;
            end
        end
        i = i + 1;
    end
end
%% plot

for iw = 1:numFiles

    for iN = 1:4

        subplot(2, 2, iN);
        plot(k, W_norm(iw).data(:, iN)); % xq(ix).data(1,1:200)-xq(ix).data(1,1)
        hold on;

        if iw == 1
            plot(k, Baselines(:, iN), 'LineWidth', 4, 'Color', 'r');
        end
        ylabel('NMF');
        ylim([0 1.1]);
    end

end

sgtitle('Subject non-negative factors obtain by averaging the experimentally obtained profiles across conditions(speed, weight)');

%% training the Speed and Weight parameters

% training the Speed parameters
v = [25, 35, 45, 55, 65, 25, 35, 45, 55, 65, 25, 35, 45, 55, 65, 25, 35, 45, 55, 65, 25, 35, 45, 55, 65];  %
for k = 1:25
    dv(k, :, :) = Increment_Speed(k).data;
end
%
[dim1, dim2, dim3] = size(dv);

%
pv_coefficients = zeros(3, dim2, dim3);

%
for i = 1:dim2
    for j = 1:dim3
        %
        dv_series = squeeze(dv(:, i, j));

        % fitting by the parameters
        pv = polyfit(v, dv_series, 2);

        % save the parameters at p
        pv_coefficients(:, i, j) = pv;
    end
end

% training the Weight parameters
w = [5, 10, 15, 5, 10, 15, 5, 10, 15, 5, 10, 15, 5, 10, 15];
% % w = [5, 7.5, 10, 12.5, 15, 5, 7.5, 10, 12.5, 15, 5, 7.5, 10, 12.5, 15, 5, 7.5, 10, 12.5, 15, 5, 7.5, 10, 12.5, 15];
selected_k = [2, 4, 7, 9, 12, 14, 17, 19, 22, 24];
index = 1;
for k = 1:25
    if ~ismember(k, selected_k)
        dw(index, :, :) = Increment_Weight(k).data; % y
        index = index + 1;
    end
    % dw(k, :, :) = Increment_Weight(k).data;
end
%
[dim1, dim2, dim3] = size(dw);

%
pw_coefficients = zeros(3, dim2, dim3);

for i = 1:dim2
    for j = 1:dim3
        %
        dw_series = squeeze(dw(:, i, j));

        % fitting by the parameters
        pw = polyfit(w, dw_series, 2);

        % save the parameters at p
        pw_coefficients(:, i, j) = pw;
    end
end

%% save the Weighting parameters to the the files

filename_pv = 'pv_coefficients.mat';
filename_pw = 'pw_coefficients.mat';

% save the data in the baseline_file
save(filename_pv, 'pv_coefficients');
save(filename_pw, 'pw_coefficients');
