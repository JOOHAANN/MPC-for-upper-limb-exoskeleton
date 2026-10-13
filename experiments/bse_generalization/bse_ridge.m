function output = bse_ridge(action, varargin)
%BSE_RIDGE Trial-balanced ridge regression with training-only scaling.
% Fit objective: mean_trial(mean_sample(error^2)) + lambda*||beta||^2.
    switch action
        case 'fit'
            x = varargin{1};
            y = varargin{2};
            group = varargin{3};
            lambda = varargin{4};
            ids = unique(group);
            weights = zeros(size(y));
            for k = 1:numel(ids)
                mask = group == ids(k);
                weights(mask) = 1 / (numel(ids) * sum(mask));
            end
            mu = weights' * x;
            scale = sqrt(weights' * (x-mu).^2);
            scale(scale < 1e-8) = 1;
            z = [ones(size(x, 1), 1), (x-mu)./scale];
            penalty = lambda * eye(size(z, 2));
            penalty(1, 1) = 0; % The intercept is not penalized.
            beta = (z' * (weights .* z) + penalty) \ (z' * (weights .* y));
            output = struct('mu', mu, 'scale', scale, 'beta', beta, 'lambda', lambda);
        case 'predict'
            model = varargin{1};
            x = varargin{2};
            output = [ones(size(x, 1), 1), (x-model.mu)./model.scale] * model.beta;
        otherwise
            error('Unknown ridge action.');
    end
end
