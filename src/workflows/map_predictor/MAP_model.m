function [weights, XP, MEP] = MAP_model(speed, elevation, baseline_mode, baseline_file)
%MAP_MODEL Predict muscle activation profiles from excitation synergies.
%   SM uses a subject-specific baseline file; GM uses the generic baseline.
%   The historical elevation argument represents load in the adapted model.
%   Original implementation: Jose Gonzalez (2014); adapted by Youhan Li.
    narginchk(3, 4);
    if nargin < 4
        if strcmp(baseline_mode, 'SM')
            error('Exo:MissingBaseline', 'SM mode requires a baseline MAT file.');
        end
        baseline_file = -1;
    end

    % Preserve this workflow's original primitive-fitting mode.
    XP = getXP_without_training(); % Reuse archived fitted primitives.
    weights = getWeights(speed, elevation, baseline_mode, baseline_file);

    % Preallocate the original 200-sample movement cycle.
    MEP = zeros(200, size(weights, 2));
    for sample = 1:200
        MEP(sample, :) = getMAP(XP, weights, sample);
    end
end
