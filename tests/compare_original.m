function compare_original(originalRoot)
%COMPARE_ORIGINAL Compare representative numerical outputs with the archive.
%   originalRoot is the original 'YOUHAN LI code and project description' folder.
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(root);
    setup_project();
    previousDir = pwd;
    cleanup = onCleanup(@() cd(previousDir));

    cd(fullfile(originalRoot, 'MAP predictor'));
    clear MAP_model getMAP getWeights getXP_without_training lowpassfilter;
    [w0, xp0, a0] = MAP_model(45, 10, 'SM', 'baseline_data.mat');
    signal = sin((0:999)' / 30) + 0.1 * cos((0:999)' / 3);
    filter0 = lowpassfilter(signal, 500, 3);
    workDir = prepare_workflow('map_predictor');
    cd(workDir);
    clear MAP_model getMAP getWeights getXP_without_training lowpassfilter;
    [w1, xp1, a1] = MAP_model(45, 10, 'SM', 'baseline_data.mat');
    assert(isequaln(w0, w1) && isequaln(xp0, xp1) && isequaln(a0, a1));
    assert(isequaln(filter0, lowpassfilter(signal, 500, 3)));

    cd(fullfile(originalRoot, 'Good demo for training BSE model', ...
        'Step 1 from Torque to determinate the ODE model'));
    clear Torque_function_constant;
    data = load('Combination.mat', 'Combination');
    samples = 1:min(500, size(data.Combination, 1));
    inputs = {zeros(1, numel(samples)), data.Combination(samples, 3), ...
        data.Combination(samples, 2), 0.001, 22.9819, 0.345, ...
        0.756835, 0.950692, 92.7424, 0.7};
    model0 = Torque_function_constant(inputs{:});
    workDir = prepare_workflow('demo_calibration__geometry');
    cd(workDir);
    clear Torque_function_constant;
    model1 = Torque_function_constant(inputs{:});
    assert(all(isfinite(model0), 'all'));
    assert(isequaln(model0, model1));
    fprintf('PASS: exact original/new equality for MAP, filtering, and 500-sample geometry dynamics.\n');
end
