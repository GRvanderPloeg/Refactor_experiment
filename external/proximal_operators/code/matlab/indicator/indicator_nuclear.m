function p = indicator_nuclear(x, eta)
%function p = indicator_nuclear(x, eta)
%
% This procedure evaluates the indicator function of the constraint set:
%
%                         ||x||_N <= eta
%
%  INPUTS
% ========
%  x   - MxN matrix
%  eta - positive, scalar or matrix compatible with the size of 'x'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (10-06-2021)
% Author  : Donnart Clara
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check the constraint
mask = (fun_nuclear(x,1) <= eta);

% evaluate the indicator function
if all(mask(:))
	p = 0;
else
	p = Inf;
end