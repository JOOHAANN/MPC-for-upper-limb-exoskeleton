% Read-only audit of the standard Demo/BSE execution path and frozen models.
folder=fileparts(mfilename('fullpath'));root=fileparts(folder);
bseFolder=fullfile(root,'Demo','BSE');addpath(bseFolder);
fprintf('Main model: %s\n',which('bse_model'));
fprintf('Main wrapper: %s\n',which('Torque_function'));
subjects=bse_subjects({'528','617','725','730'});
rows=struct([]);
for k=1:numel(subjects)
    subject=subjects{k};
    for source=1:2
        model=subject.models{source};r=model.neuralParameters;
        peak=0;driveChange=0;assistance=0;caseCount=0;
        for group=1:numel(subject.sequences)
            sequence=subject.sequences{group};
            inputs=bse_inputs(sequence,model,source);
            out=Torque_function(model,subject.balanceController,sequence,inputs,true);
            moving=out(:,20)<1e-8;
            peak=max(peak,max(abs(out(moving,22))));
            driveChange=max(driveChange,max(abs(out(moving,18:19)-inputs(moving,:)),[],'all'));
            assistance=max(assistance,max(abs(out(:,25))));
            caseCount=caseCount+numel(sequence.segments);
        end
        assert(assistance==0 && caseCount==9);
        assert(r(4)>0 && r(5)>0 && r(7)>0 && peak>0 && driveChange>0);
        names={'Experimental','MAP'};
        row=struct('Subject',string(subject.shortName),'Source',string(names{source}), ...
            'MovementKp',r(4),'MovementKd',r(5),'MovementFeedbackBlend',r(7), ...
            'Cases',caseCount,'MovementCorrectionRequestPeakNm',peak, ...
            'MaximumCommandDeviationFromInput',driveChange,'AssistancePeakNm',assistance);
        rows=[rows;row]; %#ok<AGROW>
    end
end
audit=struct2table(rows);disp(audit);
writetable(audit,fullfile(folder,'Demo_BSE_PD_audit.csv'));
fprintf('PASS: standard BSE retains active movement PD in all eight routes, across nine cases per route, with zero external assistance.\n');
