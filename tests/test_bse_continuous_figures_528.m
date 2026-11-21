function test_bse_continuous_figures_528
% Verify the continuous replay and the reported nine-case correlation means.
root=fileparts(fileparts(mfilename('fullpath')));
folder=fullfile(root,'Demo','BSE','Direct_Activation_No_Limits');
addpath(folder);output=fullfile(folder,'output');
current=load(fullfile(output,'BSE_results.mat'));
data=load(fullfile(folder,'data','subject_528.mat'));subject=data.subject;
result=current.results;

for group=1:3
    original=subject.sequences{group};
    assert(isequaln(result.subjects{1}.groups{group}.sequence,original));
    for source=1:2
        actual=result.subjects{1}.groups{group}.runs{source};
        assert(isequaln(actual.inputs,prepare_activation(original,subject.models{source},source)));
        assert(isequaln(actual.parameters(:),subject.models{source}.parameters(:)));
        assert(actual.completedSamples==numel(original.time) && all(isfinite(actual.output),'all'));
    end
end

example=result.subjects{1}.continuousExample;
q=example.sequence;parts=example.parts;
assert(numel(q.segments)==3 && all(abs(diff(q.time)-q.dt)<1e-12));
loads=cellfun(@(x)x.loadKg,q.segments);
speeds=cellfun(@(x)x.speedCode,q.segments);
assert(isequal(loads,[5 10 15]) && numel(unique(speeds))==3);
caseScore=nan(1,9);
for k=1:9
    m=result.metrics(result.metrics.CaseNumber==k,:);
    caseScore(k)=mean([m.TorqueCorrelation;m.ThetaCorrelation]);
end
cases=result.subjects{1}.cases;best=-inf;
for a=find(cellfun(@(x)x.loadKg==5,cases))
    for b=find(cellfun(@(x)x.loadKg==10,cases))
        for c=find(cellfun(@(x)x.loadKg==15,cases))
            chosen=[a b c];speedCodes=cellfun(@(x)x.speedCode,cases(chosen));
            if numel(unique(speedCodes))==3,best=max(best,mean(caseScore(chosen)));end
        end
    end
end
assert(abs(example.selectionMeanR-best)<1e-12,'The selected combination is not the best eligible one.');

% Every selected input and reference sample is copied from its original
% sequence. Joins are inside standing intervals, not inside a movement.
fields={'mass','waiting','referenceTheta','referenceTorque','torqueValid','referenceFilled'};
for k=1:numel(parts)
    p=parts{k};original=subject.sequences{p.originalGroup};
    source=(p.originalStartIndex:p.originalEndIndex)';
    target=(p.startIndex:p.endIndex)';
    for j=1:numel(fields)
        name=fields{j};assert(isequaln(q.(name)(target),original.(name)(source)));
    end
    for method=1:2
        originalInput=result.subjects{1}.groups{p.originalGroup}.runs{method}.inputs;
        assert(isequaln(example.runs{method}.inputs(target,:),originalInput(source,:)));
        assert(isequaln(example.runs{method}.parameters(:),subject.models{method}.parameters(:)));
    end
    if k<numel(parts)
        assert(q.waiting(target(end))>.5 && q.waiting(target(end)+1)>.5);
    end
end
for source=1:2
    run=example.runs{source};out=run.output;moving=q.waiting<.5;
    assert(run.completedSamples==numel(q.time) && all(isfinite(out),'all'));
    assert(all(out(:,24:26)==0,'all') && all(out(moving,20:22)==0,'all'));
    assert(max(abs(out(moving,18:19)-run.inputs(moving,:)),[],'all')<1e-12);
end

% Recalculate the displayed continuous-case R values from the raw replay.
names={'TorqueRExperimental','TorqueRMAP','ThetaRExperimental','ThetaRMAP'};
for k=1:3
    seg=q.segments{k};ix=(seg.startIndex:seg.endIndex)';
    for source=1:2
        out=example.runs{source}.output(ix,:);observed=~logical(q.referenceFilled(ix));
        valid=observed & logical(q.torqueValid(ix)) & isfinite(q.referenceTorque(ix)) & isfinite(out(:,3));
        r=corrcoef(out(valid,3),q.referenceTorque(ix(valid)));
        assert(abs(r(1,2)-example.caseMetrics.(names{source})(k))<1e-12);
        valid=observed & isfinite(q.referenceTheta(ix)) & isfinite(out(:,1));
        r=corrcoef(out(valid,1),q.referenceTheta(ix(valid)));
        assert(abs(r(1,2)-example.caseMetrics.(names{source+2})(k))<1e-12);
    end
end

summary=result.subjects{1}.correlationSummary;
for source=1:2
    m=result.metrics(result.metrics.Source==summary.Source(source),:);
    assert(height(m)==9 && summary.TorqueCases(source)==9 && summary.ThetaCases(source)==9);
    assert(abs(summary.TorqueMeanR(source)-mean(m.TorqueCorrelation,'omitnan'))<1e-12);
    assert(abs(summary.ThetaMeanR(source)-mean(m.ThetaCorrelation,'omitnan'))<1e-12);
end
for stem={'528_continuous_torque_theta','528_mean_correlation'}
    for ext={'.png','.pdf','.fig'}
        assert(isfile(fullfile(output,[stem{1} ext{1}])));
    end
end
fprintf('PASS: source replays, selection, continuous inputs/references, standing joins, fixed models and both sets of Pearson R values verified.\n');
end
