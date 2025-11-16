 function p = fun_Elliot(x)
%function p = fun_Elliot(x)
%
% This procedure evaluates the function:
%
%                  / -|x|-ln(1-|x|)-(x^2)/2  if |x| < 1
%           f(x) = | 
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
p=max(-abs(x)-log(1-abs(x))-0.5*x.^2,Inf*(abs(x)>=1));