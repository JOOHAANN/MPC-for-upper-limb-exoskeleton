function ids = list_workflows()
%LIST_WORKFLOWS List isolated experiment configurations.
    root = setup_project();
    files = dir(fullfile(root, 'config', 'workflows', '*.json'));
    ids = erase(string({files.name})', '.json');
    if nargout == 0
        fprintf('%s\n', ids);
    end
end
