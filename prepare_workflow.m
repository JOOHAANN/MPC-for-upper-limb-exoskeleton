function workDir = prepare_workflow(workflowId)
%PREPARE_WORKFLOW Create an isolated, writable experiment directory.
%   workDir = prepare_workflow('map_predictor'); cd(workDir); test_MAP
%   Each call makes a fresh directory under results. Edit canonical source
%   files, not the generated copies. Training outputs remain in this folder.
    root = setup_project();
    assert((ischar(workflowId) && isrow(workflowId)) || ...
        (isstring(workflowId) && isscalar(workflowId)), ...
        'Exo:InvalidWorkflow', 'Workflow ID must be a text scalar.');
    workflowId = char(workflowId);
    assert(~isempty(regexp(workflowId, '^[a-z0-9_]+$', 'once')), ...
        'Exo:InvalidWorkflow', 'Use an ID from list_workflows.');
    configFile = fullfile(root, 'config', 'workflows', [workflowId '.json']);
    assert(isfile(configFile), 'Exo:UnknownWorkflow', ...
        'Unknown workflow: %s. Run list_workflows.', workflowId);
    config = jsondecode(fileread(configFile));

    % Validate all inputs before making a working directory.
    for k = 1:numel(config.files)
        source = fullfile(root, config.files(k).source);
        assert(isfile(source), 'Exo:MissingInput', 'Missing input: %s', source);
    end
    parent = fullfile(root, 'results', workflowId);
    if ~isfolder(parent)
        mkdir(parent);
    end
    workDir = tempname(parent);
    mkdir(workDir);
    for k = 1:numel(config.files)
        entry = config.files(k);
        destination = fullfile(workDir, entry.target);
        if ~isfolder(fileparts(destination))
            mkdir(fileparts(destination));
        end
        copyfile(fullfile(root, entry.source), destination);
    end
    copyfile(configFile, fullfile(workDir, 'workflow_manifest.json'));
    fprintf('Prepared %s\nWorking directory: %s\n', workflowId, workDir);
    fprintf('Run cd(workDir), then the desired entry script.\n');
end
