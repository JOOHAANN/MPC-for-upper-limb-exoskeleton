% Development runner. The delivered comparison has its own readable entry.
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
addpath(fullfile(root,'Demo','BSE'),fullfile(root,'Demo','NMPC'));
out = fullfile(root,'Demo','NMPC','output','pid_nmpc_528');
if ~isfolder(out), mkdir(out); end
useNative = build_bse_solver;
subjects = bse_subjects('528'); subject = subjects{1};
settings = load(fullfile(root,'Demo','NMPC','data','controller_528.mat'));
assert(strcmp(settings.modelSHA256,subject.modelSHA256));
for group = 1:3
    sequence = subject.sequences{group};
    for source = 1:2
        name = fullfile(out,sprintf('nmpc_group%d_source%d.mat',group,source));
        if isfile(name), fprintf('Reusing %s\n',name); continue, end
        model = subject.models{source};
        plant = bse_model(model,subject.balanceController,sequence.standingAngle,useNative);
        inputs = bse_inputs(sequence,model,source);
        if source == 1, planned = bse_inputs(sequence,model,2); else, planned = inputs; end
        options = settings.controller;
        prediction = options.sourcePredictionModels{source};
        options.inputForecast = prediction.inputForecast;
        options.predictionErrorGuardNm = prediction.guardNm;
        options.predictionErrorGuardRelative = prediction.relativeGuard;
        baseline = plant.simulate(inputs,sequence.mass,sequence.waiting, ...
            sequence.initialState,sequence.dt,options.actualSubsteps);
        threshold = zeros(size(sequence.time)); limits = zeros(1,3);
        for stage = 1:3
            segment = sequence.segments{stage};
            index = segment.startIndex:segment.endIndex;
            limits(stage) = .75*max(baseline(index,4));
            threshold(index) = limits(stage);
        end
        threshold(logical(sequence.waiting)) = min(limits);
        fprintf('\nNMPC comparison cache | group %d | source %d\n',group,source);
        run = simulate_nmpc(plant,sequence,inputs,planned,threshold,options);
        modelSHA256 = subject.modelSHA256;
        save(name,'run','baseline','threshold','limits','options','modelSHA256','-v7');
        fprintf('Saved %s\n',name);
    end
end
