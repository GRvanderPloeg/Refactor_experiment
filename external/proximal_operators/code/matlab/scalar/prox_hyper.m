 function p = prox_hyper(x, gamma, d)
%function p = prox_hyper(x, gamma, d)
%
% The function computes the proximity operator of the function:
%
%              f(x) = gamma * sqrt(x^2 + d^2)
%
% When the input 'x' is an array and d, gamma are positive scalars.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  d     - positive, scalar or ND array with the same size as 'x'
% 
%  DEPENDENCIES
% ==============
%  newton.m - located in the folder 'utils'
%
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

abs_x = abs(x);

% prepare the Newton's method
fun = @(t)   t.^4 + (-2.*abs_x) .* t.^3 +   (abs_x.^2 - gamma.^2 + d.^2) .* t.^2 + (-2.*abs_x.*d.^2).*t +d.^2 .*abs_x.^2;
der = @(t) 4*t.^3 + 3*(-2.*abs_x).* t.^2 + 2*(abs_x.^2 - gamma.^2 + d.^2) .* t    + (-2.*abs_x.*d.^2);
    
% initialize the solution
p_init = abs_x./2;
    
% use the Newton's method
p = newton(fun, der, p_init);

p = p.*sign(x);