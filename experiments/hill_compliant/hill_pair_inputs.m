function inputs = hill_pair_inputs(trial, source, generator, parameters)
    % HILL_PAIR_INPUTS Experimental pair or MAP back plus generated front activation.
    % The auxiliary generator is fitted only to front-activation labels.
    if source == 1
        front = trial.front;
    else
        a = trial.a(:, 2);
        fast = memory(a, .1);
        slow = memory(a, .35);
        m = double(trial.m) / 15;
        n = double(trial.loadKg) / 15;
        if isscalar(n)
            n = ones(size(a)) * n;
        end
        assert(iscolumn(n) && numel(n) == numel(a), 'Nominal load must be scalar or a per-sample column.');
        X = [a, a.^2, fast, slow, a - slow, m, n, a .* n, a .* m, (a - fast) .* n];
        z = [ones(size(a)), (X - generator.mu(:)') ./ generator.scale(:)'] * generator.beta(:);
        front = 1 ./ (1 + exp(-max(-30, min(30, z))));
    end
    if nargin >= 4
        assert(numel(parameters) == 32 && all(isfinite(parameters)), 'Expected the complete 32-parameter Hill model.');
        back = trial.a(:, source);
        fast = memory(back, .10);
        slow = memory(back, .35);
        d = fast - slow;
        recruitment = parameters(15) * front(:) + parameters(16) * back ...
            + parameters(17) * max(0, d) + parameters(18) * max(0, -d);
        front = -expm1(-max(0, recruitment));
    end
    inputs = [trial.a(:, source), front(:)];
end

function y = memory(a, tau)
    r = exp(-.01 / tau);
    y = filter(1 - r, [1 -r], a, r * a(1));
end
