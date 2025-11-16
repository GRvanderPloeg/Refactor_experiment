 function p = prox_Logdet_NuclearNorm(x,mu,gamma)
%function p = prox_Logdet_NuclearNorm(x,mu,gamma)
%
% This procedure computes the proximity operator of the function gamma*h , and h is defined as:
%

%           h(x) = f(x) + mu*R_1(x)
% Where f(x) is defined as
%           f(x) = /   -log(Det(x))      if x is positive definite matrix
%                  \    + \infty            otherwise
% And R_1(x) the 1-Schatten norm for matrices.

%
% When the input 'x' is a square matrix n x n, the output 'p' is a vector (n).
%
%  INPUTS
% ========
%  x     - ND array ,square matrix
%  mu    - scalar
%  gamma - Real number strictly positive
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

% spectral decomposition
if size(x,2) > size(x,1)
    [u,s,v] = sv_dec_fat(x);
else
    [u,s,v] = sv_dec(x);
end

g=0.5*(s-gamma*mu+sqrt((s-mu*gamma).^2+4*gamma)) ;
% spectral reconstruction
p = sv_rec(u, g, v);
 end