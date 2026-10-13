function trials = build_bse_trials(root)
%BUILD_BSE_TRIALS Synchronize raw sensor channels without fitted MAP inputs.
% No trial-wise amplitude normalization, torque-derived activation, target
% lags, or future target samples are used as model predictors.
    addpath(root);
    setup_project();
    raw = exo_data_root();
    codes = {'5', '75', '10', '125', '15'};
    loads = [5, 7.5, 10, 12.5, 15];
    speeds = [25, 35, 45, 55, 65];
    dt = 0.02; % 50 Hz grid; interpolation does not create independent trials.
    trials = struct([]);
    for i = 1:numel(loads)
        for j = 1:numel(speeds)
            folder = sprintf('dynamics_stoop_%skg_%d', codes{i}, speeds(j));
            base = fullfile(raw, '240528Pilot', folder);
            emg = read_sto(fullfile(base, 'emgFilt.sto'), 9);
            kin = read_sto(fullfile(base, 'ik.sto'), 4);
            torque = read_sto(fullfile(base, 'Torque.sto'), 2);
            first = max([emg(1, 1), kin(1, 1), torque(1, 1)]);
            last = min([emg(end, 1), kin(end, 1), torque(end, 1)]);
            t = (0:dt:(last - first))';
            assert(numel(t) > 100, 'Insufficient common time span: %s', folder);
            e = interp1(emg(:, 1) - first, emg(:, 2:9), t, 'linear');
            q = interp1(kin(:, 1) - first, kin(:, 2:4), t, 'linear');
            y = interp1(torque(:, 1) - first, torque(:, 2), t, 'linear');
            % Fixed causal smoothers. Original acquisition filters and linear
            % timestamp interpolation still make this an offline benchmark.
            e = smooth_causal(e, dt, 0.06);
            q = smooth_causal(q, dt, 0.08);
            dq = smooth_causal([zeros(1, 3); diff(q) / dt], dt, 0.08);
            ddq = smooth_causal([zeros(1, 3); diff(dq) / dt], dt, 0.10);
            lag = max(1, round(0.08 / dt));
            ePast = [repmat(e(1, :), lag, 1); e(1:end-lag, :)];
            keep = t >= 0.5 & all(isfinite([q, dq, ddq, e, y]), 2);
            k = numel(trials) + 1;
            trials(k).id = folder;
            trials(k).loadKg = loads(i);
            trials(k).speedLabel = speeds(j);
            trials(k).time = t(keep);
            trials(k).q = q(keep, :);
            trials(k).dq = dq(keep, :);
            trials(k).ddq = ddq(keep, :);
            trials(k).emg = e(keep, :);
            trials(k).emgPast = ePast(keep, :);
            trials(k).torque = y(keep);
            trials(k).source = ['240528Pilot/' folder];
            trials(k).dt = dt;
        end
    end
end

function data = read_sto(filename, columns)
    assert(isfile(filename), 'Missing recording: %s', filename);
    fid = fopen(filename, 'r');
    cleanup = onCleanup(@() fclose(fid));
    while true
        line = fgetl(fid);
        assert(ischar(line), 'Missing endheader in %s', filename);
        if strcmp(strtrim(line), 'endheader')
            break;
        end
    end
    fgetl(fid); % Column names; do not trust archived nRows=0 headers.
    values = textscan(fid, repmat('%f', 1, columns), ...
        'Delimiter', '\t', 'MultipleDelimsAsOne', true, 'CollectOutput', true);
    data = values{1};
    data = data(all(isfinite(data), 2), :);
    [~, uniqueRows] = unique(data(:, 1), 'sorted');
    data = data(uniqueRows, :);
    assert(size(data, 2) == columns && size(data, 1) > 10);
end

function y = smooth_causal(x, dt, timeConstant)
    factor = exp(-dt / timeConstant);
    y = zeros(size(x));
    y(1, :) = x(1, :);
    for k = 2:size(x, 1)
        y(k, :) = factor * y(k-1, :) + (1-factor) * x(k, :);
    end
end
