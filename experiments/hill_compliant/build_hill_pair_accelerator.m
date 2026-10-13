function build_hill_pair_accelerator()
    % BUILD_HILL_PAIR_ACCELERATOR Optional C acceleration of the readable MATLAB solver.
    folder = fileparts(mfilename('fullpath'));
    mex('-R2017b', '-outdir', folder, '-output', 'hill_pair_mex', fullfile(folder, 'hill_pair_core.c'));
end
