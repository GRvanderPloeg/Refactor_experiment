 function p = prox_tanhA(x)
%function p = prox_tanhA(x)
%
% This procedure computes the proximity operator of the function:

%                  / x*arctan(x)+(log(1-x^2)-x^2)/2   if |x| < 1
%           f(x) = |
%                  \ +Inf                 otherwise
% When the input 'x' is an array, the output is computed element-wise.

%  INPUTS
% ========
%  x     - ND array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (10-02-2019)
% Author  : MOHAMED KERROUMI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the prox
p=tanh(x);
end
    