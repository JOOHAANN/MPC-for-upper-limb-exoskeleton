function root = setup_project()
%SETUP_PROJECT Add only project utilities, avoiding conflicting model variants.
%   Run this once after opening MATLAB in the repository root.
    root = fileparts(mfilename('fullpath'));
    addpath(root, fullfile(root, 'tools'));
end
