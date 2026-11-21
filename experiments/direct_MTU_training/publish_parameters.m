% Publish only the new direct-activation parameters, never the legacy Demo.
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');
names={'528','617','725','730'};reports=cell(4,2);parameterRows=struct([]);
parameterNames={'backCapacity','backArmToLength','backReferenceLength','jointDamping', ...
    'backMaximumVelocity','backActivationExponent','backTendonStrain','backActivationTime', ...
    'frontCapacity','frontActivationExponent','frontArmToLength','frontReferenceLength', ...
    'loadLever','loadOffset','frontInputGain','backCoactivationGain','risingGain','fallingGain', ...
    'bodyMass','bodyLever','backMomentArm','frontMomentArm','backReferenceMTLength', ...
    'frontReferenceMTLength','frontMaximumVelocity','frontTendonStrain','frontActivationTime','referencePennation'};
for k=1:numel(names)
    filename=['subject_' names{k} '.mat'];
    data=load(fullfile(folder,'before_training','data',filename));old=data.subject;
    subject=struct('id',old.id,'shortName',old.shortName,'models',{old.models}, ...
        'sequences',{old.sequences},'balanceController',old.balanceController, ...
        'activationSources',{{'Experimental activation','MAP activation'}}, ...
        'modelVersion','direct_activation_dual_MTU_no_limits');
    for source=1:2
        label=sprintf('%s_%d',names{k},source);
        reportFile=fullfile(folder,['refined_' label '.json']);
        assert(isfile(reportFile),'A complete refined fit is required: %s',label);
        report=jsondecode(fileread(reportFile));
        assert(report.dt==.01 && report.metrics.complete);
        p=report.parameters(:)';assert(numel(p)==28 && all(isfinite(p)));
        subject.models{source}.parameters=p;
        subject.models{source}.fittingWeights=report.objectiveWeights;
        reports{k,source}=report;
        for j=1:28
            row=struct('Subject',string(names{k}),'Source',source, ...
                'Parameter',string(parameterNames{j}),'Value',p(j));
            parameterRows=[parameterRows;row]; %#ok<AGROW>
        end
    end
    subject.parameterPolicy='New subject/source-specific calibration jointly fitted over nine recorded bouts; one vector per source.';
    subject.calibrationProtocol=struct('movementPD',false,'angleBoundaries',false, ...
        'assistanceNm',0,'standingPD',true,'outerStepSeconds',.01,'internalStepSeconds',.001, ...
        'independentValidation',false,'numberOfCoefficients',28,'numberOfRecordedBouts',9, ...
        'referenceUse','Training objective and metrics only; never a movement-time model input.');
    subject.calibrationProtocol.sourceWeights={reports{k,1}.objectiveWeights,reports{k,2}.objectiveWeights};
    assert(isequaln(subject.sequences,old.sequences),'Recorded inputs or references changed.');
    assert(isequaln(subject.balanceController,old.balanceController),'Standing balance was retuned.');
    save(fullfile(demo,'data',filename),'subject','-v7');
end
writetable(struct2table(parameterRows),fullfile(demo,'data','trained_parameters.csv'));
save(fullfile(demo,'data','training_reports.mat'),'reports','parameterNames','-v7');
file=fopen(fullfile(demo,'data','training_reports.json'),'w');
fprintf(file,'%s',jsonencode(reports,'PrettyPrint',true));fclose(file);
fprintf('Published eight new parameter vectors; original Demo/BSE and recorded arrays are unchanged.\n');
