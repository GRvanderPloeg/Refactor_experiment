 function p = fun_barrier3(x,c,alpha)
%function p = fun_barrier3(x,c,alpha)
%
% This procedure evaluates the function:
%
%                  /-log(alpha-||x-c||²)      if    ||x-c||²<alpha
%           f(x) = | 
%                  \ +Inf                           otherwise
%
% When the input 'x' is an array, the output 'p' is a scalar.
%
%  INPUTS
% ========
%  x     - ND array 
%  c     - ND array has the same size as x
%  alpha - scalar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (23-03-2019)
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
% check input
if  prod(size(x))~=size(c)
    error('''c'' must have the same size as ''x''');
end

p=max(-log(alpha-norm(x-c)^2),Inf*(norm(x-c)^2-alpha))  ;
 end
