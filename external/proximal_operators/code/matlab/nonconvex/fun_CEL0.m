function p = fun_CEL0(x,gamma,ai,lamb)
%function p = fun_CEL0(x,gamma,ai,lamb)
%
% This procedure evaluates the continuous exact l0 (CEL0) penalty [1].
%
% f(x) = gamma *(lamb - 0.5*ai^2*(|x| - sqrt(2*lamb)/ai)^2*1_{|x| <= sqrt(2*lamb)/ai}  
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  NOTE
% ======
% For the CEL0 penalty, ai is not a free parameter. It corresponds  to the
% l2-norm of the columns of the matrix A that is involved in the associated 
% regularized least-square problem
%  (1)    arg min 0.5||Ax-y||^2 + f(x)/gamma
% With this choice, the CEL0 relaxation (1) is an exact relaxation of the
% l2-l0 minimization problem
%  (2)    arg min 0.5||Ax-y||^2 + lamb*|x|_0
% in the sense that it preserves its global minimizers while removing some
% local minimizers [1].%
%
%  INPUTS
% ========
%  x      - ND array
%  gamma  - positive, scalar or ND array with the same size as 'x'
%  ai     - positive ND array with the same size as 'x' (SEE NOTE)
%  lamb   - positive scalar
%
%  REFERENCE
% ===========
% [1] E. Soubies, L. Blanc-Feraud, G. Aubert. A Continuous Exact l0 penalty
%     (CEL0) for Least Squares Regularized Problem.
%     SIAM Journal on Imaging Science. 8(3):1607-1639, 2015.

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
if any( ai(:) <= 0 ) || ~isscalar(ai) && any(size(ai) ~= size(x))
    error('''ai'' must be positive and of the same size as ''x''')
end
if ~isscalar(lamb) ||  lamb<=0
    error('''lamb'' must be a positive scalar')
end
%-----%

% evaluate the function
bound=sqrt(2*lamb)./ai;
ai2=ai.^2;
coef=0.5*ai2;
p =  sum(gamma(:).*(lamb - coef(:).*(abs(x(:)) - bound(:)).^2.*(abs(x(:)) <= bound(:))));

end


