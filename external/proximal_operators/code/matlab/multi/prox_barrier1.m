 function p = prox_barrier1(x,a,b,gamma)
%function p = prox_barrier1(x,a,b,gamma)
%
% This procedure computes the proximity operator of the function gamma*f , and f is defined as:

%                  /-log(b-a'*u)          if a'*u < b
%           f(x) = | 
%                  \ +Inf                 otherwise
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  INPUTS
% ========
%  x     - ND array 
%  a     - ND array has the same size as x
%  b     - Real number
%  gamma - Real number strictly positive
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
    error('''a'' must have the same size as ''x''')
end
if  (~isscalar(b))
    error('''b'' must be a scalar ''')
end
if  (gamma<=0) 
    error('''gamma''  must be strictly positive''')
end

% compute the prox
p=x+(b-a'*x-sqrt((b-a'*x)^2+4*gamma*(norm(a))^2))*a/(2*norm(a)^2);
end
    