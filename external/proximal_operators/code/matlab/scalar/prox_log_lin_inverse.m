 function p = prox_log_lin_inverse(x, gamma, k, w, r)
%function p = prox_log_lin_inverse(x, gamma, k, w, r)
%
% This procedure computes the proximity operator of the function:
%
%           / gamma * (-k*log(x) + w*x + r/x)  if x > 0
%   f(x) = |                                       with (k,w,r) >= 0 
%           \ +inf                    otherwise
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  w     - positive, scalar or ND array with the same size as 'x'
% 
%  DEPENDENCIES
% ==============
%  solver3.m - located in the folder 'utils'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (27-04-2017)
% Author  : Giovanni Chierchia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check input
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
if any( w(:) < 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
    error('''w'' must be nonnegative and either scalar or the same size as ''x''')
end
if any( r(:) < 0 ) || ~isscalar(r) && any(size(r) ~= size(x))
    error('''r'' must be nonnegative and either scalar or the same size as ''x''')
end
if any( k(:) < 0 ) || ~isscalar(k) && any(size(k) ~= size(x))
    error('''k'' must be nonnegative and either scalar or the same size as ''x''')
end
%-----%


% compute the prox
p =  solver3(1, gamma.*w-x, -gamma.*k, -r.*gamma);