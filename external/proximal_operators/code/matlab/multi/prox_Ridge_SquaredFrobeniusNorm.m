 function p = prox_Ridge_SquaredFrobeniusNorm(x,mu,gamma)
%function p = prox_Ridge_SquaredFrobeniusNorm(x,mu,gamma)
%
% This procedure computes the proximity operator of the function gamma*f , and f is defined as:
%

%           f(x) = (1/2)*||x||_F^2 + mu*||x||_F^2

%
% When the input 'x' is a square matrix n x n, the output 'p' is a vector (n).
%
%  INPUTS
% ========
%  x     - ND array ,square matrix
%  mu    - scalar
%  gamma - Real number strictly positive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (30-03-2019)
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

g=s/(1+gamma*(1+2*mu))  ;

% spectral reconstruction
p = sv_rec(u, g, v);
 end