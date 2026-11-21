function inputs = prepare_activation(sequence, model, source)
%PREPARE_ACTIVATION Use recorded abdomen or its independently trained MAP.
% MAP scaling was fixed from experimental mean curves, separately per load.
back = sequence.backActivation(:, source);
if source == 1
    abdomen = sequence.abdomenActivation;
else
    assert(isfield(sequence,'abdomenMAPActivation'), ...
        'The abdominal MAP must be trained and packaged first.');
    abdomen = sequence.abdomenMAPActivation;
end
assert(numel(model.parameters)==24,'This demo requires the strict 24-parameter model.');
assert(model.abdomenScale==1,'Additional fitted input gains are disabled.');
inputs = max(0,min(1,[back abdomen]));
end
