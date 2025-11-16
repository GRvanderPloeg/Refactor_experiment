function p = fun_SCAD(x,gamma,a,lamb)
%function p = fun_SCAD(x,gamma,a,lamb)
%
% This procedure evaluates the smoothly clipped absolute deviation (SCAD) 
% penalty. [1]
%
%                |  lamb*|x|                                if |x| <= lamb
% f(x) = gamma * |  -(lamb^2-2*a*lamb*|x| + x^2)/[2*(a-1)]  if lamb < |x| <= a*lamb   
%                |  (a+1)*lamb^2/2                          if |x| > a*lamb
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  INPUTS
% ========
%  x      - ND array
%  gamma  - positive, scalar or ND array with the same size as 'x'
%  a      - > 2, scalar or ND array with the same size as 'x'
%  lamb   - positive, scalar or ND array with the same size as 'x'
%
%  REFERENCE
% ===========
% [1] J. Fan and R. Li, Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties.
%     Journal of the American Statistical Association, 96(456):1348-1360, 2001.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (21-10-2019)
% Author  : Emmanuel Soubies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
if any( a(:) <= 2 ) || ~isscalar(a) && any(size(a) ~= size(x))
    error('''a'' must be lerger than 2 and either scalar or the same size as ''x''')
end
if any( lamb(:) <= 0 ) || ~isscalar(lamb) && any(size(lamb) ~= size(x))
    error('''lamb'' must be positive and either scalar or the same size as ''x''')
end
%-----%

% evaluate the function
p = sum(gamma(:).*(   lamb(:).*abs(x(:))                                            .* (abs(x(:)) <=lamb(:)) ...
                - (lamb(:).^2 - 2*a(:).*lamb(:).*abs(x(:)) + x(:).^2)./(2*(a(:)-1)) .* (abs(x(:)) <= a(:).*lamb(:)).*(abs(x(:)) >lamb(:)) ...
                + (a(:)+1).*lamb(:).^2/2                                            .* (abs(x(:)) > a(:).*lamb(:))));

end

