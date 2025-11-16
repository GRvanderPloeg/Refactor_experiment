 function p = fun_UniSig(x)
%function p = fun_UniSig(x)
%
% This procedure evaluates the function:
%
%                  /(x+1/2)*ln(x+1/2)+(1/2 -x)*log(1/2 -x)-(x^2+1/4)/2  if |x| < 1/2
%           f(x) = | -1/4                 |x| = 1/2
%                  \ +Inf                 otherwise
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
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
%-----%
p=max(min((abs(x)<0.5).*((x+0.5).*log(x+0.5)+(0.5-x).*log(0.5-x)-0.5*(x.^2+0.25)),-0.25*(abs(x)==0.5)),Inf*(abs(x)>0.5));