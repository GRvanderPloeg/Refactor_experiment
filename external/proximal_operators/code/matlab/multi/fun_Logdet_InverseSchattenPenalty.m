 function p = fun_Logdet_InverseSchattenPenalty(x,mu,p)
%function p = fun_Logdet_InverseSchattenPenalty(x,mu,p)
%
% This procedure evaluates the function:
%

%           h(x) = f(x) + mu*R_p(x^-1)^p
% Where f(x) is defined as
%           f(x) = /   -log(Det(x))      if x is positive definite matrix
%                  \    + \infty             otherwise
% And R_p(x^-1) the p-Schatten norm of the inverse of x.
% When the input 'x' is a square matrix n x n, the output 'p' is a scalar.
%
%  INPUTS
% ========
%  x     - ND array ,square matrix
%  mu    - scalar
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
if  p<1
    error('''p'' must be >= 1');
end
if  det(x)==0
    error('''x'' must be invertible');
end
y=inv(x);
s=eig(y);
p= -log(det(x))+mu*sum(abs(s).^p)   ;
 end
