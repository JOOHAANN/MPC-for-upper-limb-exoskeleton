function inputs = bse_inputs(sequence, model, source)
%BSE_INPUTS Prepare experimental or trained-MAP excitation for the Hill pair.
% backActivation(:,2) contains the frozen trained MAP output, not raw EMG.
% The abdominal MAP input uses the fitted activation-only regression.
back = sequence.backActivation(:, source);
if source == 1
    abdomen = sequence.abdomenActivation;
else
    fast = memory(back, .10, sequence.dt);
    slow = memory(back, .35, sequence.dt);
    load = sequence.nominalLoad / 15;
    mass = sequence.mass / 15;
    features = [back, back.^2, fast, slow, back-slow, mass, load, ...
        back.*load, back.*mass, (back-fast).*load];
    generator = model.generator;
    design = [ones(size(back)), ...
        (features-generator.mu(:)')./generator.scale(:)'];
    logit = design * generator.beta(:);
    abdomen = 1 ./ (1 + exp(-max(-30, min(30, logit))));
end
p = model.parameters;
% Gains 15:18 recruit abdomen from its input, back input and activation changes.
change = memory(back, .10, sequence.dt) - memory(back, .35, sequence.dt);
recruitment = p(15)*abdomen + p(16)*back ...
    + p(17)*max(change, 0) + p(18)*max(-change, 0);
inputs = [back, -expm1(-max(recruitment, 0))];
end

function value = memory(signal, tau, dt)
decay = exp(-dt/tau);
value = filter(1-decay, [1 -decay], signal, decay*signal(1));
end
