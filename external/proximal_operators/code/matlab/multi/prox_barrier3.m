 function p = prox_barrier3(x,c,alpha,gamma)
%function p = prox_barrier3(x,c,alpha,gamma)
%
% This procedure computes the proximity operator of the function gamma*f , and f is defined as:

%                  /-log(alpha-||x-c||²)      if    ||x-c||²<alpha
%           f(x) = | 
%                  \ +Inf                           otherwise
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  INPUTS
% ========
%  x     - ND array 
%  c     - ND array has the same size as x
%  alpha - scalar
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
if  prod(size(x))~=size(c)
    error('''c'' must have the same size as ''x''');
end
if  (gamma<=0) 
    error('''gamma''  must be strictly positive''');
end
% compute the prox
A= roots([1,-norm(x-c),-(alpha+2*gamma),alpha*norm(x-c)]);
A = A(imag(A)==0);
p=c+(alpha-A(0<=A & A<sqrt(alpha))^2)*(x-c)/(alpha- A(0<=A & A<sqrt(alpha))^2+2*gamma);
end
    