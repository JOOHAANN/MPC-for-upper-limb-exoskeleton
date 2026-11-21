% Publish the recalibrated 528 vectors; keep experimental load scaling fixed.
support=fileparts(mfilename('fullpath'));root=fileparts(fileparts(support));
folder=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');
data=load(fullfile(support,'calibration_input.mat'));subject=data.subject;
selection=jsondecode(fileread(fullfile(support,'selected_fits.json')));
reports=cell(1,2);rows=struct([]);
names={'backCapacity','backArmToLength','backReferenceLength','jointDamping', ...
    'backMaximumVelocity','backActivationExponent','backTendonStrain','backActivationTime', ...
    'abdomenCapacity','abdomenActivationExponent','abdomenArmToLength','abdomenReferenceLength', ...
    'loadLever','loadOffset','bodyMass','bodyLever','backMomentArm','abdomenMomentArm', ...
    'backReferenceMTLength','abdomenReferenceMTLength','abdomenMaximumVelocity', ...
    'abdomenTendonStrain','abdomenActivationTime','referencePennation'};
for source=1:2
    report=jsondecode(fileread(fullfile(support,selection.files{source})));
    assert(report.dt==.01 && report.metrics.complete);
    p=report.parameters(:)';assert(numel(p)==24 && report.abdomenScale==1);
    subject.models{source}.parameters=p;subject.models{source}.abdomenScale=1;
    reports{source}=report;
    for k=1:24
        row=struct('Subject',"528",'Source',source,'Parameter',string(names{k}),'Value',p(k));
        rows=[rows;row]; %#ok<AGROW>
    end
end
subject.calibrationProtocol.status='calibrated';
subject.calibrationProtocol.outerStepSeconds=.01;subject.calibrationProtocol.internalStepSeconds=.001;
subject.calibrationProtocol.numberOfInputGains=0;
subject.calibrationProtocol.referenceUse='Training objective and evaluation only; never movement feedback.';
save(fullfile(folder,'data','subject_528.mat'),'subject','-v7');
writetable(struct2table(rows),fullfile(folder,'data','trained_parameters.csv'));
save(fullfile(folder,'data','training_reports.mat'),'reports','-v7');
file=fopen(fullfile(folder,'data','training_reports.json'),'w');
fprintf(file,'%s',jsonencode(reports,'PrettyPrint',true));fclose(file);
fprintf('Published two fixed 24-parameter BSE vectors. MAP amplitude calibration remains fixed by load.\n');
