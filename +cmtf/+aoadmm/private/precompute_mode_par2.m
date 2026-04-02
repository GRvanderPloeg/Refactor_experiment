function [A_m, B_m, C_m, rho_m, L_m, G, inner_iters, last_mttkrp, last_had, last_m] = ...
    precompute_mode_par2(Z, G, G_transp_G, m, p, iter, coupl_id, ...
                         last_mttkrp, last_had, last_m, options)
% Precompute A, B, C, rho (and Cholesky L) for one PAR2 mode.
%
% Branches on the mode position within the PARAFAC2 model:
%   Mode 1 (coupled mode)  — matrix A/B/C/rho; L_m = []; G unchanged.
%   Mode 2 ("funny" mode)  — per-slice cell A/B/C/rho/L; G unchanged.
%   Mode 3 (diagonal mode) — per-slice cell A/B/C/rho; for unconstrained
%                            uncoupled modes the row-wise ALS update of
%                            G.fac{m} is done here (inner_iters = 1);
%                            for constrained uncoupled modes L_m{k} is
%                            built so the caller can invoke dispatch_uncoupled.
%
% inner_iters is [] unless this function performs the actual factor update
% (mode 3 unconstrained uncoupled), in which case it returns 1.

    mode_pos = find(Z.modes{p} == m);
    L_m = [];
    inner_iters = [];

    switch mode_pos
        case 1  % first PARAFAC2 mode (can be coupled)
            A_m = zeros(size(G.fac{m}));
            C_m = zeros(size(G.fac{m}, 2), size(G.fac{m}, 2));
            for k = 1:length(Z.size{Z.modes{p}(2)})
                A_m = A_m + Z.object{p}{k} * G.fac{Z.modes{p}(2)}{k} * diag(G.fac{Z.modes{p}(3)}(k,:));
                C_m = C_m + diag(G.fac{Z.modes{p}(3)}(k,:)) * G_transp_G{Z.modes{p}(2)}{k} * diag(G.fac{Z.modes{p}(3)}(k,:));
            end
            last_had{p} = C_m;
            last_mttkrp{p} = A_m;
            last_m(p) = 1;
            A_m = Z.weights(p) * A_m;
            rho_m = trace(C_m) / size(C_m, 1);
            B_m = Z.weights(p) * C_m;
            if isfield(Z, 'ridge')
                B_m = B_m + Z.ridge(m) * eye(size(B_m));
            end
            if options.bsum
                A_m = A_m + options.bsum_weight/2 * G.fac{m};
                B_m = B_m + options.bsum_weight/2 * eye(size(B_m));
            end

        case 2  % second PARAFAC2 mode (the "funny" mode, no external coupling)
            K = length(Z.size{Z.modes{p}(2)});
            A_m = cell(K, 1);
            B_m = cell(K, 1);
            C_m = cell(K, 1);
            rho_m = zeros(K, 1);
            L_m = cell(K, 1);
            for k = 1:K
                A_m{k} = Z.weights(p) * Z.object{p}{k}' * G.fac{Z.modes{p}(1)} * diag(G.fac{Z.modes{p}(3)}(k,:));
                C_m{k} = diag(G.fac{Z.modes{p}(3)}(k,:)) * G_transp_G{Z.modes{p}(1)} * diag(G.fac{Z.modes{p}(3)}(k,:));
                rho_m(k) = trace(C_m{k}) / size(C_m{k}, 1);
                if isfield(options, 'increase_factor_rhoBk')
                    rho_m(k) = options.increase_factor_rhoBk * rho_m(k);
                end
                B_m{k} = Z.weights(p) * C_m{k};
                B_m{k} = B_m{k} + rho_m(k)/2 * eye(size(B_m{k}));  % always coupled
                if isfield(Z, 'ridge')
                    B_m{k} = B_m{k} + Z.ridge(m) * eye(size(B_m{k}));
                end
                if options.bsum
                    A_m{k} = A_m{k} + options.bsum_weight/2 * G.fac{m}{k};
                    B_m{k} = B_m{k} + options.bsum_weight/2 * eye(size(B_m{k}));
                end
                if Z.constrained_modes(m) && iter >= options.iter_start_PAR2Bkconstraint
                    B_m{k} = B_m{k} + rho_m(k)/2 * eye(size(B_m{k}));
                end
                L_m{k} = chol(B_m{k}, 'lower');
            end
            last_m(p) = 2;

        otherwise  % third PARAFAC2 mode (diagonal, can be coupled)
            K = length(Z.size{Z.modes{p}(2)});
            A_m = cell(K, 1);
            B_m = cell(K, 1);
            C_m = cell(K, 1);
            rho_m = zeros(K, 1);
            L_m = cell(K, 1);
            for k = 1:K
                A_m{k} = Z.weights(p) * diag(G.fac{Z.modes{p}(1)}' * Z.object{p}{k} * G.fac{Z.modes{p}(2)}{k});
                C_m{k} = G_transp_G{Z.modes{p}(1)} .* G_transp_G{Z.modes{p}(2)}{k};
                rho_m(k) = trace(C_m{k}) / size(C_m{k}, 1);
                B_m{k} = Z.weights(p) * C_m{k};
                if isfield(Z, 'ridge')
                    B_m{k} = B_m{k} + Z.ridge(m) * eye(size(B_m{k}));
                end
                %last_mttkrp{p}{k} = A_m{k} * 1/Z.weights(p);
                if options.bsum
                    A_m{k} = A_m{k} + options.bsum_weight/2 * G.fac{m}(k,:)';
                    B_m{k} = B_m{k} + options.bsum_weight/2 * eye(size(B_m{k}));
                end
                if coupl_id == 0
                    if Z.constrained_modes(m) == 0  % unconstrained: do ALS update here
                        G.fac{m}(k,:) = (B_m{k} \ A_m{k})';
                    else  % constrained: build Cholesky for dispatch_uncoupled
                        B_m{k} = B_m{k} + rho_m(k)/2 * eye(size(B_m{k}));
                        L_m{k} = chol(B_m{k}', 'lower');
                    end
                end
            end
            last_m(p) = 3;
            if coupl_id == 0 && Z.constrained_modes(m) == 0
                inner_iters = 1;
            end
    end
end
