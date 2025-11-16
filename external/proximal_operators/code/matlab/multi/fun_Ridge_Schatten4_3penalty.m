 function p = fun_Ridge_Schatten4_3penalty(x,mu)
%function p = fun_Ridge_Schatten4_3penalty(x,mu)
%
% This procedure evaluates the function:
%

%           f(x) = (1/2)*||x||² + mu*R_4/3^(4/3)(x)
%           with R_p(x) the p-Shatten norm for matrices.

%
% When the input 'x' is a square matrix n x n, the output 'p' is a scalar.
%
%  INPUTS
% ========
%  x     - ND array ,square matrix
%  mu    - scalar
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
s=eig(x);
p=0.5*norm(x,'fro')^2+mu*sum(abs(s).^(4/3))  ;
 end
