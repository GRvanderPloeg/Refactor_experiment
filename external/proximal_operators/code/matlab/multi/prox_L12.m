function [sol] = prox_l12(x, gamma, W )
%function p = prox_l12(x, gamma)
%
% This procedure computes the proximity operator of the function
%
%                    f(x) = gamma * ||x||_{12}
%
% where || x ||_12 =  sqrt ( sum_j ( sum_i |x(i,j)|)^2  )
%
%  INPUTS
% ========
%  x     - matrix
%  gamma - positive, scalar (or ND array compatible with the size of 'x')
%  DEPENDENCIES
% ==============
%  prox_L1.m - located in the folder 'multi'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (11-06-2021)
% Author  : Donnart Clara/ Nathanael Perraudin (epfl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Handle the parameter p
% argmin_x 1/2 ||x-y||_2^2 + p/gamma || x ||_{qp}^p
gamma = 2*gamma;

z = transpose(x);

if nargin<3
    w = ones(size(z));
else
    w = transpose(W);
end

r = abs(z./w);

[Nelg, Ng] = size(r);

[r, ind] = sort(r,1,'descend');

zo = zeros(Nelg,Ng);
wo = zeros(Nelg,Ng);

for jj = 1 : Ng
    zo(:,jj) = z(ind(:,jj),jj);
    wo(:,jj) = w(ind(:,jj),jj);
end

Mg = zeros(Ng,1);
Kw = zeros(Ng,1);
ny = zeros(Ng,1);

for jj = 1 : Ng

    for ii=1:Nelg-1
        temp = gamma * ...
            sum( wo(1:(ii+1),jj).^2 .* ( r(1:(ii+1),jj) - r((ii+1),jj) ) )...
            - r((ii+1),jj);
        if temp >= 0
            Mg(jj) = ii;
            Kw(jj) = sum( wo(1:ii,jj).^2 );
            ny(jj) = norm( zo(1:ii,jj),1);
            break;
        end
    end
    
    % handle exception
    if Mg(jj)==0
        Mg(jj) = 1;
        Kw(jj) = wo(1,jj).^2;
        ny(jj) = abs(zo(1,jj));
    end

end

tau = gamma./( 1 + repmat(Kw', size(x,2),1) * gamma) .* repmat(ny', size(x,2),1);
sol = prox_L1( X, tau)';
end



