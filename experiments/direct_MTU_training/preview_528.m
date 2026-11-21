% Render the current fitted 528 candidate at the full MATLAB replay step.
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');
data=load(fullfile(folder,'before_training','data','subject_528.mat'));subject=data.subject;
for source=1:2
    file=fullfile(folder,sprintf('refined_528_%d.json',source));
    if ~isfile(file),file=fullfile(folder,sprintf('refined_528_%d_0.02.json',source));end
    report=jsondecode(fileread(file));subject.models{source}.parameters=report.parameters(:)';
end
subject.parameterPolicy='New subject/source-specific direct-MTU calibration; candidate under full-step verification.';
subject.calibrationProtocol='Movement PD and angle boundaries absent; standing balance unchanged.';
save(fullfile(demo,'data','subject_528.mat'),'subject','-v7');
Subject_list={'528'};Show_figures=false;Use_native=true;Output_label='528_candidate';
run(fullfile(demo,'run_dual_MTU_BSE.m'));
