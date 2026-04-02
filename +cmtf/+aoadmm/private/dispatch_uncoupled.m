function [G, G_transp_G, sum_column_norms_sqr, inner_iters, lbfgsb_iters] = ...
    dispatch_uncoupled(Z, G, nb_modes, G_transp_G, sum_column_norms_sqr, ...
                       A_m, B_m, L_m, rho_m, m, p, rho, lscalar, uscalar, fh, gh, options)
% Handle the coupl_id==0 (uncoupled) factor update for mode m.
%
% Performs the unconstrained ALS update or constrained ADMM update for an
% uncoupled mode, then updates G_transp_G{m} (Frobenius) or
% sum_column_norms_sqr(m) (non-Frobenius).
%
% Inputs:
%   A_m, B_m  - precomputed MTTKRP and Gram matrix for this mode
%   L_m       - Cholesky factor if already built (pass [] to build here)
%   rho_m     - scalar ADMM penalty for this mode
%   rho       - full rho cell (passed to ADMM_constrained_only)
%   All other inputs are standard problem/solver fields.

    lbfgsb_iters = [];

    if Z.constrained_modes(m) == 0  % unconstrained
        if strcmp(Z.loss_function{p}, 'Frobenius')
            G.fac{m} = A_m / B_m;
        else
            [lbfgsb_iters, G] = lbfgsb_update(Z, G, lscalar, uscalar, fh, gh, p, m, false, -1, rho_m, options);
        end
        inner_iters = 1;
    else  % constrained — use ADMM
        if strcmp(Z.loss_function{p}, 'Frobenius') && isempty(L_m)
            % Build Cholesky on a local copy; does not modify caller's B_m
            B_augmented = B_m + rho_m/2 * eye(size(B_m));
            L_m = chol(B_augmented', 'lower');
        end
        [inner_iters, lbfgsb_iters, G] = ADMM_constrained_only(Z, G, nb_modes, lscalar, uscalar, fh, gh, A_m, L_m, m, p, rho, options);
    end

    % Update cached Gram matrix or column norms after factor update
    if strcmp(Z.loss_function{p}, 'Frobenius')
        G_transp_G{m} = G.fac{m}' * G.fac{m};
    else
        for r = 1:size(G.fac{m}, 2)
            sum_column_norms_sqr(m, 1) = norm(G.fac{m}(:, r))^2;
        end
    end
end
