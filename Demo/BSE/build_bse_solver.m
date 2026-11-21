function available = build_bse_solver(force)
%BUILD_BSE_SOLVER Compile the optional shared Hill accelerator once.
% The readable MATLAB equations in bse_model remain available as a fallback.
if nargin < 1
    force = false;
end
folder = fileparts(mfilename('fullpath'));
source = fullfile(folder, 'private', 'native', 'bse_hill_mex.c');
target = fullfile(folder, 'private', ['bse_hill_mex_v2.' mexext]);
sources = dir(fullfile(folder, 'private', 'native', '*'));
newest = max([sources(~[sources.isdir]).datenum]);
binary = dir(target);
if ~force && ~isempty(binary) && binary.datenum >= newest
    available = true;
    return
end
if isempty(mex.getCompilerConfigurations('C', 'Selected'))
    available = false;
    warning('BSE:Compiler', ['No C compiler is configured. The MATLAB solver is available, ' ...
        'but NMPC will be slow. Configure a compiler with mex -setup C.']);
    return
end
fprintf('Building the shared Hill solver...\n');
clear bse_hill_mex_v2
mex('-R2017b', '-O', '-outdir', fullfile(folder, 'private'), ...
    '-output', 'bse_hill_mex_v2', source);
available = true;
end
