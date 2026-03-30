function [Z] = project_ortho(X)
% Project a matrix onto the set of matrices with orthonormal columns via SVD.
%
% Computes the nearest matrix with orthonormal columns to X by computing
% the economy SVD X = U*S*V' and returning Z = U*V'.
%
% Syntax:
%   Z = cmtf.project.project_ortho(X)
%
% Input:
%   X - Matrix of size (m x n), m >= n
%
% Output:
%   Z - Matrix of size (m x n) with orthonormal columns (Z'*Z = I)
[U,~,V] = svd(X,0);
Z = U*V';
end

