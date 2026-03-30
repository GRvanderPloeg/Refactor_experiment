function setup()
%SETUP  Add the CMTF Toolbox and its external dependencies to the MATLAB path.
%
% Adds the toolbox root and any present external libraries (Tensor Toolbox,
% Poblano, L-BFGS-B wrapper, proximal operators) to the path. Silently skips
% dependencies whose folders are not found.
%
% Syntax:
%   setup()
%
% TODO: detect whether compile_mex.m needs to be run.

root = fileparts(mfilename('fullpath'));

% Add main code
addpath(genpath(fullfile(root, 'src'))); % or +cmtf if already packaged
addpath(fullfile(root));  % top-level functions, if any

% Add external dependencies
external = fullfile(root, 'external');

% Tensor Toolbox
if isfolder(fullfile(external, 'tensor_toolbox'))
    addpath(genpath(fullfile(external, 'tensor_toolbox')));
end

% Poblano
if isfolder(fullfile(external, 'poblano_toolbox'))
    addpath(genpath(fullfile(external, 'poblano_toolbox')));
end

% LBFGSB
if isfolder(fullfile(external, 'lbfgsb_wrapper'))
    addpath(genpath(fullfile(external, 'lbfgsb_wrapper')));
end

% Proximal operators
if isfolder(fullfile(external, 'proximal_operators'))
    addpath(genpath(fullfile(external, 'proximal_operators')));
end

% Optional: additional helper libs
% addpath(genpath(fullfile(external,'something_else')));

fprintf('[CMTF] Setup complete. Dependencies added to path.\n');
end
