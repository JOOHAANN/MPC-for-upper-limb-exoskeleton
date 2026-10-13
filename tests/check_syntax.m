function report = check_syntax()
%CHECK_SYNTAX Collect MATLAB Code Analyzer diagnostics without running models.
    root = fileparts(fileparts(mfilename('fullpath')));
    files = [dir(fullfile(root, 'src', '**', '*.m')); ...
             dir(fullfile(root, 'tools', '*.m')); ...
             dir(fullfile(root, '*.m')); dir(fullfile(root, 'tests', '*.m'))];
    report = struct('file', {}, 'diagnostics', {});
    for k = 1:numel(files)
        filename = fullfile(files(k).folder, files(k).name);
        report(k).file = erase(filename, [root filesep]);
        report(k).diagnostics = checkcode(filename, '-id');
    end
    out = fullfile(root, 'results', 'code_analyzer.json');
    fid = fopen(out, 'w', 'n', 'UTF-8');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', jsonencode(report, PrettyPrint=true));
    fprintf('Analyzed %d MATLAB files. Report: %s\n', numel(files), out);
end
