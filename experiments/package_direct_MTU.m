% Package existing recordings and physical coefficients without neural models.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');
for name={'528','617','725','730'}
    file=['subject_' name{1} '.mat'];
    data=load(fullfile(root,'Demo','BSE','data',file));subject=data.subject;
    for source=1:2
        old=subject.models{source};
        subject.models{source}=struct('parameters',old.parameters(1:28),'generator',old.generator);
    end
    subject.parameterPolicy='Transferred physical/recruitment coefficients; not refitted after movement-feedback removal.';
    subject.modelVersion='direct_activation_dual_MTU_no_angle_limits';
    save(fullfile(folder,'data',file),'subject','-v7');
end
fprintf('Packaged four subjects. Neural/intent and angle-boundary coefficients are absent from the new model vectors.\n');
