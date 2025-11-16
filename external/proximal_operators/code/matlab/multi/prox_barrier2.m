 function p = prox_barrier2(x,a,bm,bM,gamma)
%function p = prox_barrier2(x,a,bm,bM,gamma)
%
% This procedure computes the proximity operator of the function gamma*f , and f is defined as:

%                  /-log(bM-a'*u)-log(a'*u-bm)      if    bm<a'*u < bM
%           f(x) = | 
%                  \ +Inf                           otherwise
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
    error('''a'' must have the same size as ''x''');
end
if  (~isscalar(bm)) | (~isscalar(bM))
    error('''bm'' and ''bM'' must be a scalar ''');
end
if  bm>=bM
    error('''you must have bm<bM''');
end
if  (gamma<=0) 
    error('''gamma''  must be strictly positive''');
end
% compute the prox
A= roots([1,-(bm+bM+a'*x),bm*bM+(a'*x)*(bm+bM)-2*gamma*(norm(a)^2),-bm*bM*a'*x+gamma*(bm+bM)*norm(a)^2]);
A = A(imag(A)==0);
p=x+(A(bm<A & A<bM)-a'*x)*a/(norm(a)^2);
end
    