% Package true experimental abdominal input and its independently trained MAP.
folder=fileparts(mfilename('fullpath'));root=fileparts(fileparts(folder));
demo=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');addpath(demo);
data=load(fullfile(folder,'mixed_abdomen_candidate','data','subject_528.mat'));
subject=data.subject;map=load(fullfile(demo,'data','abdominal_MAP_528.mat'));
edges=[];
for k=1:numel(map.records)
    r=map.records{k};profile=r.profile(:);
    edges=[edges;profile(1:5);profile(end-4:end)]; %#ok<AGROW>
end
baseline=median(edges);
for source=1:2
    old=subject.models{source}.parameters;
    subject.models{source}=struct('parameters',old([1:14 19:28]), ...
        'abdomenScale',1,'inputPolicy','Recorded abdomen or independently trained abdominal MAP; constant amplitude scaling only.');
end
for group=1:numel(subject.sequences)
    q=subject.sequences{group};q.abdomenMAPActivation=baseline*ones(size(q.time));
    for stage=1:numel(q.segments)
        seg=q.segments{stage};ix=seg.startIndex:seg.endIndex;
        duration=q.time(ix(end))-q.time(ix(1));
        phase=(q.time(ix)-q.time(ix(1)))/duration;
        q.abdomenMAPActivation(ix)=predict_abdominal_MAP(map.model,seg.loadKg,60/duration,phase);
    end
    subject.sequences{group}=q;
end
subject.modelVersion='strict_recorded_abdomen_dual_MTU_no_limits';
subject.parameterPolicy='One fixed 24-parameter physical vector and one abdominal amplitude gain per activation source, jointly fitted over nine bouts.';
subject.abdominalMAP=map.model;
subject.calibrationProtocol=struct('movementPD',false,'standingPD',true, ...
    'angleBoundaries',false,'assistanceNm',0,'numberOfPhysicalCoefficients',24, ...
    'numberOfInputGains',1,'numberOfRecordedBouts',9,'independentValidation',false);
save(fullfile(demo,'data','subject_528.mat'),'subject','-v7');
save(fullfile(folder,'strict_input_528.mat'),'subject','-v7');
fprintf('Packaged strict recorded/MAP abdominal inputs for subject 528.\n');
