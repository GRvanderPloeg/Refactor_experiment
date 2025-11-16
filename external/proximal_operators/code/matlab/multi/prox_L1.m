function p = prox_L1(x, gamma)
%function p = prox_L1(x, gamma)
%
% This procedure computes the proximity operator of the function
%
%                    f(x) = gamma * ||x||_1
%
%  INPUTS
% ========
%  x     - ND array or MxN matrix
%  gamma - positive, scalar (or ND array/matrix compatible with the size of 'x')
%  DEPENDENCIES
% ==============
%  prox_abs.m - located in the folder 'multi'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (02-06-2021)
% Author  : Donnart Clara
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = size(x);

if any( gamma(:) <= 0 ) || ~isscalar(gamma) && (any(size(gamma)~=sz))
    error('''gamma'' must be positive and either scalar or compatible with the blocks of ''x''')
end
%-----%

p = prox_abs(x, gamma);

