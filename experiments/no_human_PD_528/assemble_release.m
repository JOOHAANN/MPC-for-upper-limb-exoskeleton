% Assemble independently simulated continuous sequences without changing data.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
output=fullfile(root,'Demo','NMPC','output');
destination=fullfile(output,'no_human_PD_528');
for group=1:3
    sourceFolder=fullfile(folder,'development_runs',sprintf('no_human_PD_release_group%d',group));
    for source=1:2
        filename=sprintf('group%d_source%d.mat',group,source);
        cache=load(fullfile(sourceFolder,filename));
        assert(numel(cache.sequence.time)>4000,'A shortened probe is not a release.');
        copyfile(fullfile(sourceFolder,filename),fullfile(destination,filename));
    end
end
