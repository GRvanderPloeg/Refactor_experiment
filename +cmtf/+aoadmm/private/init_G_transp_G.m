function [G_transp_G, sum_column_norms_sqr] = init_G_transp_G(Z, G, which_p, nb_modes)
% Initialise G'*G cache and column-norm cache before the outer ADMM loop.
%
% For each mode m:
%   Frobenius CP      — G_transp_G{m} = G.fac{m}' * G.fac{m}
%   Frobenius PAR2 1  — G_transp_G{m} = G.fac{m}' * G.fac{m}
%   Frobenius PAR2 2  — G_transp_G{m}{k} = G.fac{m}{k}' * G.fac{m}{k}  (per slice)
%   Non-Frobenius     — sum_column_norms_sqr(m) accumulated from column norms

    G_transp_G = cell(nb_modes, 1);
    sum_column_norms_sqr = zeros(nb_modes, 1);

    for m = 1:nb_modes
        p = which_p(m);
        if strcmp(Z.loss_function{p}, 'Frobenius')
            if strcmp(Z.model{p}, 'CP')
                G_transp_G{m} = G.fac{m}' * G.fac{m};
            elseif strcmp(Z.model{p}, 'PAR2')
                if 1 == find(Z.modes{p} == m)
                    G_transp_G{m} = G.fac{m}' * G.fac{m};
                elseif 2 == find(Z.modes{p} == m)
                    for k = 1:length(Z.size{m})
                        G_transp_G{m}{k} = G.fac{m}{k}' * G.fac{m}{k};
                    end
                end
            end
        else
            for r = 1:size(G.fac{m}, 2)
                sum_column_norms_sqr(m, 1) = sum_column_norms_sqr(m, 1) + norm(G.fac{m}(:, r))^2;
            end
        end
    end
end
