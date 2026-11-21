% Probe the effort/motion trade-off without modifying the original controller.
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
useNative = build_bse_solver;
subjects = bse_subjects('528'); subject = subjects{1};
selection = jsondecode(fileread(fullfile(folder,'pid_selection.json')));
g = selection.gains;
pid = struct('Kp',g(1),'Ki',g(2),'Kd',g(3),'derivativeTime',g(4));
outputFolder = fullfile(root,'Demo','NMPC','output','pid_nmpc_528');
for source = 1:2
    group = 1; sequence = subject.sequences{group};
    cache = load(fullfile(outputFolder,sprintf('nmpc_group%d_source%d.mat',group,source)));
    model = subject.models{source};
    plant = bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
    inputs = bse_inputs(sequence,model,source);
    if source == 1, planned = bse_inputs(sequence,model,2); else, planned = inputs; end
    pidRun = simulate_threshold_pid(plant,sequence,inputs,cache.threshold,cache.options,pid);
    fprintf('PID source %d: effort %.4f, angle %.6f\n',source, ...
        sum(pidRun.actual(:,25).^2)*sequence.dt, ...
        sqrt(mean((pidRun.actual(:,1)-cache.baseline(:,1)).^2)));
    for weight = [0 1 100 1000]
        options = cache.options;
        options.torqueConstraintPriority = false;
        options.backTorqueMarginNm = 0;
        options.backTorqueUnderWeight = 0;
        options.assistanceEffortWeight = weight;
        options.recoveryBackTargetNm = 0;
        fprintf('\nPROBE source %d effort weight %g\n',source,weight);
        run = simulate_nmpc(plant,sequence,inputs,planned,cache.threshold,options);
        fprintf('RESULT effort %.4f angle %.6f excess %.4f\n', ...
            sum(run.actual(:,25).^2)*sequence.dt, ...
            sqrt(mean((run.actual(:,1)-cache.baseline(:,1)).^2)), ...
            max(run.actual(:,4)-cache.threshold));
        save(fullfile(folder,sprintf('probe_source%d_weight%g.mat',source,weight)), ...
            'run','options','pidRun','pid','-v7');
    end
end
