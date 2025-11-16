function p = prox_linf1(x, gamma)
%function p = prox_linf1(x, gamma)
%
% This procedure computes the proximity operator of the function
%
%                    f(x) =  gamma * \sum_{i=1}^N \sup_{1\le j \le M} |X_{i,j}|
%
%  INPUTS
% ========
%  x     - MxN matrix
%  gamma - positive, scalar (or ND array compatible with the size of 'x')
%  DEPENDENCIES
% ==============
%  prox_Linf.m - located in the folder 'multi'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (11-06-2021)
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
%-----%

%p = prox_Linf(x, gamma);
p = zeros(size(x));
if isscalar(gamma) == 1
    for i = 1:size(x,1)
        p(i,:) = prox_Linf(x(i,:),gamma,1);
    end
else
    for i = 1:size(x,1)
       p(i,:) = prox_Linf(x(i,:),gamma(i,:), 1);
    end
end



