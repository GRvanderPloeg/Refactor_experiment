 function p = fun_barrier2(x,a,bm,bM)
%function p = fun_barrier2(x,a,bm,bM)
%
% This procedure evaluates the function:
%
%                  /-log(bM-a'*u)-log(a'*u-bm)      if    bm<a'*u < bM
%           f(x) = | 
%                  \ +Inf                           otherwise
%
% When the input 'x' is an array, the output 'p' is a scalar.
%
%  INPUTS
% ========
%  x     - ND array 
%  a     - ND array has the same size as x
%  bm    - scalar
%  bM    - scalar such as bm < bM
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

p=max(max(-log(bM-a'*x)-log(a'*x-bm),Inf*(bm-a'*x)),Inf*(a'*x-bM))  ;
 end
