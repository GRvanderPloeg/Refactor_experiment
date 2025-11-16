 function p = fun_Ridge_FrobeniusNorm(x,mu)
%function p = fun_Ridge_FrobeniusNorm(x,mu)
%
% This procedure evaluates the function:
%

%           f(x) = (1/2)*||x||_F^2 + mu*||x||_F

%
% When the input 'X' is a square matrix n x n, the output 'p' is a scalar.
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
p=0.5*norm(x,'fro')^2+mu*norm(x,'fro')  ;
 end
