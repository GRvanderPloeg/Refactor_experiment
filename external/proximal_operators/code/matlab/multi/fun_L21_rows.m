function p = fun_L21_rows_Clara(x, gamma)
%function p = fun_L21_Clara(x, gamma)
%
% This procedure evaluates the function:
%
%                    f(x) = gamma * \sum_{i=1}^N|\sum_{j=1}^M|X(i,j)|^2|^{\frac{1}{2}}
%
%  INPUTS
% ========
%  x     - MxN matrix
%  gamma - positive, scalar or ND array compatible with the size of 'x'
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

for i = 1:size(x,1)
  if ~isscalar(gamma) ==0
    p = p+ fun_L2(x(i,:), gamma);
  else
    p= p+fun_L2(x(i,:), gamma(i,:));
  end
end
