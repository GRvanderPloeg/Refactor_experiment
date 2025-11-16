 function p = project_nuclear(x, eta)
%function p = project_nuclear(x, eta)
%
% This procedure computes the projection onto the constraint set:
%
%                         ||x||_N <= eta
%  INPUTS
% ========
%  x   - MxN matrix
%  eta - positive, scalar or matrix compatible with the size of 'x'
% 
%  DEPENDENCIES
% ==============
%  project_L1.m - located in the folder 'indicator'
%  prox_svd.m     - located in the folder 'utils'
%  prox_Linf.m - located in the folder 'multi'
%  prox_max.m  - located in the folder 'multi'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (10-06-2021)
% Author  : Donnart Clara
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = prox_svd(x, eta, @project_L1); 