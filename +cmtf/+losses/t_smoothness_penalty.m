function loss = t_smoothness_penalty(x, smoothness_l)
% Compute the temporal smoothness penalty for a sequence of factor matrices.
%
% Returns smoothness_l times the sum of squared Frobenius norms of differences
% between consecutive factor matrices, encouraging smooth variation across slices.
%
% Syntax:
%   loss = cmtf.losses.t_smoothness_penalty(x, smoothness_l)
%
% Inputs:
%   x            - Cell array of K factor matrices {B_1, B_2, ..., B_K}
%   smoothness_l - Non-negative scalar scaling factor for the penalty
%
% Output:
%   loss - Scalar penalty: smoothness_l * sum_{i=2}^{K} ||B_i - B_{i-1}||_F^2

    loss = 0;

    for i=2:length(x)
        loss = loss + norm(x{i} - x{i-1},'fro')^2;
    end

    loss = smoothness_l * loss;

end