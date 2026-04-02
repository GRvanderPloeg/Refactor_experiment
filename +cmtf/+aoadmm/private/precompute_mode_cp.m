function [A_m, B_m, C_m, rho_m, last_mttkrp, last_had, last_m] = ...
    precompute_mode_cp(Z, G, G_transp_G, sum_column_norms_sqr, ...
                       m, p, last_mttkrp, last_had, last_m, options)
% Precompute A, B, C, rho for one CP mode (Frobenius or non-Frobenius).
%
% For Frobenius loss: computes MTTKRP (A_m), Hadamard product of Gram
% matrices (C_m), Gram matrix (B_m), penalty parameter (rho_m), and
% updates last_mttkrp/last_had/last_m for efficient function evaluation.
% For non-Frobenius loss: only rho_m is meaningful; A_m/B_m/C_m are [].

    if strcmp(Z.loss_function{p}, 'Frobenius')
        if length(size(Z.object{p})) >= 3  % Tensor
            A_m = Z.weights(p) * mttkrp(Z.object{p}, G.fac(Z.modes{p}), find(Z.modes{p}==m));
            C_m = ones(size(G_transp_G{m}));
            for j = Z.modes{p}
                if j ~= m
                    C_m = C_m .* G_transp_G{j};
                end
            end
        else  % Matrix
            matrix_mode = find(Z.modes{p} == m);
            if matrix_mode == 1
                A_m = Z.weights(p) * double(Z.object{p}) * G.fac{Z.modes{p}(2)};
                C_m = G_transp_G{Z.modes{p}(2)};
            else
                A_m = Z.weights(p) * double(Z.object{p})' * G.fac{Z.modes{p}(1)};
                C_m = G_transp_G{Z.modes{p}(1)};
            end
        end
        rho_m = trace(C_m) / size(C_m, 1);
        B_m = Z.weights(p) * C_m;
        if isfield(Z, 'ridge')
            B_m = B_m + Z.ridge(m) * eye(size(B_m));
        end

        last_mttkrp{p} = A_m * 1/Z.weights(p);
        last_had{p} = C_m;
        last_m(p) = m;

        if options.bsum
            A_m = A_m + options.bsum_weight/2 * G.fac{m};
            B_m = B_m + options.bsum_weight/2 * eye(size(B_m));
        end
    else  % non-Frobenius
        A_m = [];
        B_m = [];
        C_m = [];
        rho_m = sum(sum_column_norms_sqr([1:m-1, m+1:end]));
    end
end
