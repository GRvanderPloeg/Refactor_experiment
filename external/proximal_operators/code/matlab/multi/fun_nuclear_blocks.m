function p = fun_nuclear_blocks(X, gamma,ind_r, ind_c, W)
%function p = fun_nuclear_blocks(x, gamma,ind_r, ind_c, W)
%
% This procedure evaluates blockwise the function:
%
%                    f(x) = gamma * \sum_{i= 1}^N \sum_{j=1}^M w_{i,j} \|X_{i,j}\|_N
% 
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array compatible with the size of 'x'
%  ind_r - Vector partitioning the rows of X in groups EXAMPLE: ind_r [1 2 2 3 3 1] 
%  means that the first block contains the first and last rows of x
%  ind_c - Vector partitioning the columns of X in groups 
%
%  DEPENDENCIES
% ==============
%  fun_svd.m         - located in the folder 'utils'
%  fun_abs.m        - located in the folder 'scalar'
%  fun_nuclear.m     - located in the folder 'multi'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (04-06-2021)
% Author  : Donnart Clara
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(X);

% check input
if any( W(:) <= 0 ) ||  any(size(W) ~= [max(ind_r), max(ind_c)])
    error('''W'' must be positive and of size (max(ind_r),max(ind_x))')
end

if length(ind_r) ~= m || length(ind_c) ~= n
    error('ind_r or ind_c has not the right length');
end


% evaluates the function
p = 0;

for i = 1: max(ind_r)
  for j = 1:max(ind_c)
    I = (ind_r == i);
    J = (ind_c == j);
    M = mtimes(I',J);
    A = X.*M;       %The elements that don't belong to the block are replaced by zeros
    p = p+ fun_nuclear(A, W(i,j)*gamma); %Update of p
  end
end