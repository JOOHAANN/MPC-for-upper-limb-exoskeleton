% Inspect optional OpenMP compilation without changing the working binary.
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
source=fullfile(root,'Demo','NMPC','private','native','bse_no_posture_mex.c');
flags={'COMPFLAGS=$COMPFLAGS /openmp','COMPFLAGS="$COMPFLAGS /openmp"'};
for k=1:numel(flags)
    fprintf('\nCompiler flag form %d\n',k);
    try
        mex('-v',flags{k},'-R2017b','-O','-outdir',folder,'-output','omp_probe',source);
        fprintf('OPENMP BUILD SUCCEEDED\n');
        break
    catch problem
        disp(getReport(problem,'extended','hyperlinks','off'));
    end
end
