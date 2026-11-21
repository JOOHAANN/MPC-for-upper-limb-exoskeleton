function probe_weights
% Joint-objective calibration on one initial action, not a held-out test.
weights=[1000 100 5000 100;500 100 20000 100;200 300 50000 100];
for k=1:size(weights,1)
    run_one(weights(k,:),k);
end
end
function run_one(weights,index)
NMPC_weights=weights; Run_label=sprintf('no_human_PD_probe%d',index);
Source_list=1; Group_list=1; Show_figures=false; Reuse_cache=false; Max_samples=1400;
run('D:/matlab/Journal Paper of Upper-Limb Exoskeleton/GitHub_Ready/Demo/NMPC/offline_PID_NMPC_no_human_PD_528.m');
end
