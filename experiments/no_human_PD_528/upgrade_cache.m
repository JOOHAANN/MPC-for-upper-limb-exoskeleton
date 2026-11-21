% One-time cache migration after PID-only calibration. Preserve provenance.
folder=fileparts(mfilename('fullpath')); root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','NMPC'); addpath(demo);
old=no_posture_nmpc_signature(demo,fullfile(folder,'initial_main.m'));
current=no_posture_nmpc_signature(demo);
assert(strcmp(old,current),'NMPC dynamics/optimization code changed; rerun instead.');
for group=1:3
    for source=1:2
        name=fullfile(demo,'output','no_human_PD_528',sprintf('group%d_source%d.mat',group,source));
        cache=load(name); cache.nmpcSignature=current;
        cache.nmpcSourceSHA256=cache.signature;
        save(name,'-struct','cache','-v7');
    end
end
