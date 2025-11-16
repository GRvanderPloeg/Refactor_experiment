 function p = prox_UniSig(x)
%function p = prox_UniSig(x)
%
% This procedure computes the proximity operator of the function:

%                  /(x+1/2)*ln(x+1/2)+(1/2 -x)*log(1/2 -x)-(x^2+1/4)/2  if |x| < 1/2
%           f(x) = | -1/4                 |x| = 1/2
%                  \ +Inf                 otherwise
% When the input 'x' is an array, the output is computed element-wise.

%  INPUTS
% ========
%  x     - ND array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (11-02-2019)
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
p=1./(1+exp(-x))-0.5;
end
    