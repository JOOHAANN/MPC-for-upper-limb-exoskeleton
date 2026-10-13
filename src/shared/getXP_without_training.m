function XP = getXP_without_training()
%GETXP_WITHOUT_TRAINING Load previously fitted excitation primitives.
%   Returns the samples-by-synergies matrix without refitting or plotting.
%   Original implementation: Jose Gonzalez (2014).
    fitted = load('Primitives_fit_y.mat', 'Primitives_fit_y');
    XP = fitted.Primitives_fit_y;
end
