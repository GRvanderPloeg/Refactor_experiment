function p = prox_MCP(x,gamma,bet,lamb)
%function p = prox_MCP(x,gamma,bet,lamb)
%
% This procedure computes the proximity operator of the minimax concave 
% penalty (MCP) [1]
%
%              | lamb^2*bet/2              if  |x| > bet*lamb
% f(x) = gamma*|                                                 
%	           | [|x| - x^2/(2*bet*lamb)]  otherwise  
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x      - ND array
%  gamma  - positive, scalar or ND array with the same size as 'x'
%  bet    - positive, scalar or ND array with the same size as 'x'
%  lamb   - positive, scalar or ND array with the same size as 'x'
%
%  REFERENCE
% ===========
% [1] C.H. Zhang. Nearly umbiased variable selection under minimax concave penalty.
%     The Annals of Statistics 38(2):894-942, 2010

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
if any( lamb(:) <= 0 ) || ~isscalar(lamb) && any(size(lamb) ~= size(x))
    error('''lamb'' must be positive and either scalar or the same size as ''x''')
end
if any( bet(:) <= 0 ) || ~isscalar(bet) && any(size(bet) ~= size(x))
    error('''bet'' must be positive and either scalar or the same size as ''x''')
end
%-----%

p =   x .* (abs(x) > sqrt(gamma.*bet.*lamb.^2))                             .* (bet <= gamma) ...
    + sign(x) .* min(bet./(bet-gamma) .* max(abs(x)-lamb.*gamma,0)  , abs(x)) .* (bet >gamma);

end


