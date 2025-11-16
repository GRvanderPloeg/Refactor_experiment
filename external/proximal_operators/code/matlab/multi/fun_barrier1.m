 function p = fun_barrier1(x,a,b)
%function p = fun_barrier1(x,a,b)
%
% This procedure evaluates the function:
%
%                  /-log(b-a'*x)          if a'*x < b
%           f(x) = | 
%                  \ +Inf                 otherwise
%
% When the input 'x' is an array, the output 'p' is a scalar.
%
%  INPUTS
% ========
%  x     - ND array 
%  a     - ND array has the same size as x
%  b     - Real number
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
if  prod(size(x))~=size(a)
    error('''a'' must have the same size as ''x''');
end
if  (~isscalar(b))
    error('''b'' must be a scalar ''');
end
p=max(-log(b-a'*x),Inf*(a'*x-b)) ;
 end
