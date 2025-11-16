 function p = fun_Indicator1 (x)
%function p = fun_Indicator1(x)
%
% This procedure evaluates the function
%
%                       f(x) = 0 if  -1<=x<=1
%                              Inf if not
% When the input 'x' is an array, the output is computed element-wise.

%  INPUTS
% ========
%  x     - ND array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (01-02-2019)
% Author  : MOHAMED KERROUMI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluate the function
p=Inf*ones(size(x));
p=p.*(x<-1 | x>1);
p(isnan(p))=0;
 end
 