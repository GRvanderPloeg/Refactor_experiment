 function p = fun_hyper(x, gamma, d)
%function p = fun_hyper(x, gamma, d)
%
% The function evaluate the function:
%
%       f(x) = gamma * sqrt(x^2 + d^2)
%
% When the input 'x' is an array and d, gamma are positive scalars.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  d     - positive, scalar or ND array with the same size as 'x'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (27-04-2017)
% Author  : Emilie Chouzenoux
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
if any( d(:) <= 0 ) || ~isscalar(d) && any(size(d) ~= size(x))
    error('''d'' must be positive and either scalar or the same size as ''x''')
end
%-----%


% evaluate the function
p = sum(gamma(:).*sqrt(x(:).^2 + d(:).^2)); 
