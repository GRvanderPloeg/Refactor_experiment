 function p = fun_PReLu(x, alpha)
%function p = fun_PReLu(x, alpha)
%
% This procedure evaluates the function:
%
%                  / 0   if x > 0
%           f(x) = |
%                  \ (1/alpha -1)*(x^2)/2  otherwise
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  INPUTS
% ========
%  x     - ND array
%  alpha - scalar in ]0,1] or ND array with the same size as 'x'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (05-02-2019)
% Author  : MOHAMED KERROUMI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% check input
if ~isscalar(alpha) && any(size(alpha) ~= size(x))
    error('''alpha'' must be either scalar or the same size as ''x''')
end
if ~all( 0< alpha(:) & alpha(:) <= 1 )
    error('''alpha'' must be in '']0,1]''')
end
%-----%
p = zeros(size(x));

% evaluate the function
p = max(0,-sign(x).*(1/alpha -1).*(x.^2)/2);