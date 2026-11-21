% Publish only the selected strict-input 528 calibration; preserve recordings.
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');
data=load(fullfile(folder,'strict_input_528.mat'));subject=data.subject;
selected=jsondecode(fileread(fullfile(folder,'strict_selection.json')));
reports=cell(1,2);rows=struct([]);
names={'backCapacity','backArmToLength','backReferenceLength','jointDamping', ...
    'backMaximumVelocity','backActivationExponent','backTendonStrain','backActivationTime', ...
    'abdomenCapacity','abdomenActivationExponent','abdomenArmToLength','abdomenReferenceLength', ...
    'loadLever','loadOffset','bodyMass','bodyLever','backMomentArm','abdomenMomentArm', ...
    'backReferenceMTLength','abdomenReferenceMTLength','abdomenMaximumVelocity', ...
    'abdomenTendonStrain','abdomenActivationTime','referencePennation','abdomenInputScale'};
for source=1:2
    report=jsondecode(fileread(fullfile(folder,selected.files{source})));
    assert(report.dt==.01 && report.metrics.complete);
    p=report.parameters(:)';assert(numel(p)==24);
    subject.models{source}.parameters=p;
    subject.models{source}.abdomenScale=report.abdomenScale;
    reports{source}=report;values=[p report.abdomenScale];
    for k=1:25
        row=struct('Subject',"528",'Source',source,'Parameter',string(names{k}),'Value',values(k));
        rows=[rows;row]; %#ok<AGROW>
    end
end
subject.calibrationProtocol.outerStepSeconds=.01;
subject.calibrationProtocol.internalStepSeconds=.001;
subject.calibrationProtocol.referenceUse='Training loss and evaluation only; never runtime movement feedback.';
subject.calibrationProtocol.abdominalInput='Recorded rect_abd_r or separately trained abdominal MAP; one amplitude gain per source.';
save(fullfile(demo,'data','subject_528.mat'),'subject','-v7');
writetable(struct2table(rows),fullfile(demo,'data','trained_parameters.csv'));
save(fullfile(demo,'data','training_reports.mat'),'reports','-v7');
file=fopen(fullfile(demo,'data','training_reports.json'),'w');
fprintf(file,'%s',jsonencode(reports,'PrettyPrint',true));fclose(file);
fprintf('Published two strict-input parameter vectors for all nine 528 conditions.\n');
