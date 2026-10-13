function build_hill_accelerator()
%BUILD_HILL_ACCELERATOR Optional acceleration; requires a supported C compiler.
    folder=fileparts(mfilename('fullpath'));
    mex('-R2017b','-outdir',folder,fullfile(folder,'hill_integrate_mex.c'));
end
