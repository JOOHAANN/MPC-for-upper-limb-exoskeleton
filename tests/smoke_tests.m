function smoke_tests()
%SMOKE_TESTS Check manifest isolation, recording paths, and MAP reconstruction.
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(root);
    setup_project();
    previousDir = pwd;
    cleanup = onCleanup(@() cd(previousDir));
    ids = list_workflows();
    assert(numel(ids) == 38, 'Unexpected workflow count.');
    for k = 1:numel(ids)
        workDir = prepare_workflow(ids(k));
        manifest = jsondecode(fileread(fullfile(workDir, 'workflow_manifest.json')));
        for j = 1:numel(manifest.files)
            assert(isfile(fullfile(workDir, manifest.files(j).target)));
        end
    end
    workDir = prepare_workflow('map_predictor');
    cd(workDir);
    [weights, XP, activation] = MAP_model(45, 10, 'SM', 'baseline_data.mat');
    assert(isequal(size(activation), [200, 6]));
    assert(all(isfinite(activation), 'all'));
    assert(max(abs(activation - XP(1:200, :) * weights), [], 'all') < 1e-12);
    % Verify that the resolver supports the legacy CalibrationData level.
    raw = exo_data_root();
    if isfolder(fullfile(raw, '240528Pilot'))
        a = exo_data_file(raw, '\240528Pilot\CalibrationData\dynamics_stoop_5kg_45\', 'emgFilt.sto');
        b = exo_data_file(raw, '240528Pilot/dynamics_stoop_5kg_45', 'emgFilt.sto');
        assert(strcmp(a, b));
    end
    fprintf('PASS: 38 workflow manifests, MAP reconstruction, and raw-data paths.\n');
end
