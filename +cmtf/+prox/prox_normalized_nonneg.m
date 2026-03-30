function [Y] = prox_normalized_nonneg(X)
% Proximal operator projecting each column onto the nonnegative unit sphere.
%
% For each column of X: zeroes out negative entries, then normalizes to unit
% L2-norm. If a column is entirely non-positive, the entry with the largest
% original value is set to 1 and the rest to 0.
%
% Syntax:
%   Y = cmtf.prox.prox_normalized_nonneg(X)
%
% Input:
%   X - Matrix of size (m x R)
%
% Output:
%   Y - Matrix of size (m x R) with each column nonneg and unit L2-norm
    Y = project_box(X,0,inf); % non-negativity
    for r=1:size(Y,2)
        if norm(Y(:,r),2)==0
            [~,maxcoord] = max(X(:,r));
            Y(maxcoord,r) = 1; %set maximum coordinate to 1, leave rest 0
        else
            Y(:,r) = Y(:,r)./norm(Y(:,r),2); %normalize
        end
    end
end

