% Calibrate one torque penalty per activation source to match total PID effort.
% Only the measured effort enters this search; angle/torque outcomes do not.
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
useNative = build_bse_solver;
subjects = bse_subjects('528'); subject = subjects{1};
selection = jsondecode(fileread(fullfile(folder,'pid_selection.json')));
g = selection.gains;
pid = struct('Kp',g(1),'Ki',g(2),'Kd',g(3),'derivativeTime',g(4));
settings = load(fullfile(root,'Demo','NMPC','data','controller_528.mat'));
outputFolder = fullfile(root,'Demo','NMPC','output','pid_nmpc_528');
cacheFolder = fullfile(folder,'tuning_cache');
if ~isfolder(cacheFolder), mkdir(cacheFolder); end
tolerance = .01;
torqueWeights = zeros(1,2);
matching = cell(1,2);
for source = 1:2
    sequences = subject.sequences;
    plants = cell(1,3); observed = cell(1,3); planned = cell(1,3);
    pidRuns = cell(1,3); thresholds = cell(1,3); baselines = cell(1,3);
    target = 0;
    for group = 1:3
        sequence = sequences{group}; model = subject.models{source};
        plants{group} = bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
        observed{group} = bse_inputs(sequence,model,source);
        if source==1, planned{group}=bse_inputs(sequence,model,2); else, planned{group}=observed{group}; end
        cache = load(fullfile(outputFolder,sprintf('nmpc_group%d_source%d.mat',group,source)));
        thresholds{group} = cache.threshold; baselines{group} = cache.baseline;
        pidRuns{group} = simulate_threshold_pid(plants{group},sequence,observed{group}, ...
            thresholds{group},settings.controller,pid);
        target = target+sum(pidRuns{group}.actual(:,25).^2)*sequence.dt;
    end
    history = [];
    lowerWeight = 0; upperWeight = inf;
    candidate = 64; bestError = inf;
    for iteration = 1:12
        options = comparison_options_528(settings.controller,source,candidate);
        runs = cell(1,3); effort = 0;
        for group = 1:3
            name = fullfile(cacheFolder,sprintf('matched_trial_source%d_weight%.9g_group%d.mat',source,candidate,group));
            if isfile(name)
                cache = load(name); run = cache.run;
            else
                fprintf('\nMATCH source %d | iteration %d | weight %.9g | group %d\n',source,iteration,candidate,group);
                run = simulate_nmpc(plants{group},sequences{group},observed{group}, ...
                    planned{group},thresholds{group},options);
                save(name,'run','options','-v7');
            end
            runs{group} = run;
            effort = effort+sum(run.actual(:,25).^2)*sequences{group}.dt;
        end
        relativeError = effort/target-1;
        fprintf('MATCH RESULT source %d weight %.9g PID %.6f NMPC %.6f relative %.6f\n', ...
            source,candidate,target,effort,relativeError);
        history = [history;candidate,target,effort,relativeError]; %#ok<AGROW>
        if abs(relativeError)<bestError
            bestError = abs(relativeError); bestRuns = runs; bestOptions = options; bestWeight = candidate;
        end
        if abs(relativeError)<=tolerance, break, end
        if relativeError<0
            lowerWeight = candidate;
            if isinf(upperWeight), candidate = 2*candidate; else, candidate = sqrt(lowerWeight*upperWeight); end
        else
            upperWeight = candidate;
            if lowerWeight==0, candidate = candidate/2; else, candidate = sqrt(lowerWeight*upperWeight); end
        end
    end
    torqueWeights(source) = bestWeight;
    matching{source} = struct('targetEffort',target,'relativeError',bestError,'history',history);
    for group = 1:3
        run = bestRuns{group}; pidRun = pidRuns{group};
        threshold = thresholds{group}; baseline = baselines{group}; options = bestOptions;
        modelSHA256 = subject.modelSHA256;
        save(fullfile(outputFolder,sprintf('comparison_group%d_source%d.mat',group,source)), ...
            'run','pidRun','threshold','baseline','options','pid','modelSHA256','-v7');
    end
    save(fullfile(folder,'matching_progress.mat'),'pid','torqueWeights','matching','-v7');
end
protocol = struct('effortDefinition','integral of squared delivered torque', ...
    'effortMatchScope','sum over all three continuous sequences (nine conditions) per source', ...
    'effortMatchTolerance',tolerance,'pidTuningGroups',1,'pidHeldOutGroups',[2 3], ...
    'nmpcMatchingUsesOnlyEffort',true,'sameTracesForRQ1AndRQ2',true);
modelSHA256 = subject.modelSHA256;
save(fullfile(root,'Demo','NMPC','data','comparison_settings_528.mat'), ...
    'pid','torqueWeights','protocol','matching','modelSHA256','-v7');
fprintf('\nFIXED PID: Kp %.12g Ki %.12g Kd %.12g Tf %.12g\n',pid.Kp,pid.Ki,pid.Kd,pid.derivativeTime);
fprintf('FIXED NMPC torque weights: %.12g %.12g\n',torqueWeights);
