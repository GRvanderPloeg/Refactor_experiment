function p = fun_L12(x, gamma)
%function p = fun_L12(x, gamma)
%
% This procedure evaluates the function:
%
%                    f(x) = gamma * ||x||_12
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array compatible with the blocks of 'x'
%
%  DEPENDENCIES
% ==============
%  fun_L2.m - located in the folder 'multi'
%
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

sz = size(x);
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && (any(size(gamma)~=sz))
    error('''gamma'' must be positive and either scalar or compatible with the blocks of ''x''')
end

p = 0;

for j = 1:size(x,1)
  if ~isscalar(gamma) ==0
    p = p+ fun_L2(x(:,j), gamma)^2;
  else
    p= p+fun_L2(x(:,j), gamma(:,j))^2;
  end
end
p = sqrt(p);