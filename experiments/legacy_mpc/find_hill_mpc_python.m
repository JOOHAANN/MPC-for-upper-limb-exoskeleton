function executable = find_hill_mpc_python()
% Locate the Python environment used by the Hill/MPC simulation and plots.
executable = getenv('EXO_PYTHON');
if ~isempty(executable), return; end
try
    settings = pyenv;
    executable = char(settings.Executable);
catch
    executable = '';
end
if ~isempty(executable), return; end
profileDirectory = getenv('USERPROFILE');
candidates = {fullfile(profileDirectory, 'miniconda3', 'python.exe'), ...
              fullfile(profileDirectory, 'anaconda3', 'python.exe')};
for index = 1:numel(candidates)
    if isfile(candidates{index})
        executable = candidates{index};
        return;
    end
end
if ispc, executable = 'python'; else, executable = 'python3'; end
end
