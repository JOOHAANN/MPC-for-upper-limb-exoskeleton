function MAP = getMAP(XP, weights, time)
%GETMAP Reconstruct muscle activation at one sample from synergy weights.
%   XP is samples-by-synergies; weights is synergies-by-muscles.
%   Original implementation: Jose Gonzalez (2014).
    MAP = XP(time, :) * weights;
end
