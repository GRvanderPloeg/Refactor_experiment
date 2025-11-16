 function p = fun_Logdet_SquaredFrobeniusNorm(x,mu)
%function p = fun_Logdet_SquaredFrobeniusNorm(x,mu)
%
% This procedure evaluates the function:
%

%           h(x) = f(x) + mu*||x||_F²
% Where f(x) is defined as
%           f(x) = /   -log(Det(x))      if x is positive definite matrix
%                  \    + \infty             otherwise
% And ||x||_F the Frobenius norm for matrices.
% When the input 'x' is a square matrix n x n, the output 'p' is a scalar.
%
%  INPUTS
% ========
%  x     - ND array ,square matrix
%  mu    - scalar
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
p= -log(det(x))+mu*norm(x,'fro')^2  ;
 end
