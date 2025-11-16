 function p = prox_Logdet_SchattenPenalty(x,mu,gamma,p)
%function p = prox_Logdet_SchattenPenalty(x,mu,gamma,p)
%
% This procedure computes the proximity operator of the function gamma*h , and h is defined as:
%

%           h(x) = f(x) + mu*R_p(x)^p
% Where f(x) is defined as
%           f(x) = /   -log(Det(x))      if x is positive definite matrix
%                  \    + \infty             otherwise
% And R_p(x) the p-Schatten norm for matrices.
% When the input 'x' is a square matrix n x n, the output 'p' is a vector (n).
%
%  INPUTS
% ========
%  x     - ND array ,square matrix
%  mu    - scalar
%  gamma - Real number strictly positive
%  p     - Non-zero integer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (02-04-2019)
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
if  mu<=0
    error('''mu'' must be stricly positive');
end
if  (gamma<=0) 
    error('''gamma''  must be strictly positive''');
end
if  p<1
    error('''p'' must be >= 1');
% spectral decomposition
if size(x,2) > size(x,1)
        [u,s,v] = sv_dec_fat(x);
else
    [u,s,v] = sv_dec(x);
end
g=zeros(1,length(s));
for i=1 :length(s)
    a=zeros(1,p+1);
    a([1,p-1,p,p+1])=[gamma*mu*p,1,-s(i),-gamma];
    r=roots(a);
    r = r(imag(r)==0);
    g(i)=r(r>0);
    
end

% spectral reconstruction
p = sv_rec(u, g, v);
 end