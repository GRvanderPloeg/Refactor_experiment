    function [f_tensors,f_couplings,f_constraints,f_PAR2_couplings] = CMTF_AOADMM_func_eval(Z, G, G_transp_G, Znorm_const,last_mttkrp,last_had,last_m, fh)
    % evaluates the 'residuals' of the objective function
    % ftensors = sum_i w_i ||T_i-[|C_i,1,C_i,2,C_i,3|]||_F^2 + g_i,d(C_i,d) (for all regularizations g_i,d)
    % f_couplings = sum_i ||C_i-Delta_i||_F^2
    % f_constraints = sum_i ||C_i-Z_i||_F^2

    P = numel(Z.object);
    nb_modes = numel(Z.size);

    fp = zeros(P,1);
    for pp = 1:P
        if strcmp(Z.model{pp},'CP')
            if strcmp(Z.loss_function{pp},'Frobenius')
                if length(size(Z.object{pp}))>=3
                    % Tensor 
                    if isempty(last_mttkrp)
                       fp(pp) = cmtf.losses.cp_func(Z.object{pp}, G.fac(Z.modes{pp}),  Znorm_const{pp},Z.weights(pp)); 
                    else
                        f_1 =  Znorm_const{pp};
                        V = last_mttkrp{pp}.*G.fac{last_m(pp)};
                        f_2 = sum(V(:));
                        W = last_had{pp}.*G_transp_G{last_m(pp)};
                        f_3 = sum(W(:));
                        f = f_1 - 2* f_2 + f_3;
                        fp(pp) = Z.weights(pp) *f;
                    end
                elseif length(size(Z.object{pp}))==2
                    % Matrix   
                   if isempty(last_mttkrp)
                        fp(pp) = cmtf.losses.pca_func(Z.object{pp}, G.fac(Z.modes{pp}),  Znorm_const{pp},Z.weights(pp));   
                    else
                        f_1 =  Znorm_const{pp};
                        V = last_mttkrp{pp}.*G.fac{last_m(pp)};
                        f_2 = sum(V(:));
                        W = last_had{pp}.*G_transp_G{last_m(pp)};
                        f_3 = sum(W(:));
                        f = f_1 - 2* f_2 + f_3;
                        fp(pp) = Z.weights(pp) *f;
                    end
                end
            else
                fp(pp) = Z.weights(pp)*(Znorm_const{pp} + collapse(fh{pp}(Z.object{pp},full(ktensor(G.fac(Z.modes{pp}))))));   % can we avoid this???????????????    
            end
        elseif strcmp(Z.model{pp},'PAR2')
            fp(pp) = 0;
            if ~isempty(last_mttkrp) && last_m(pp)==1
                f_1 =  Znorm_const{pp};
                V = last_mttkrp{pp}.*G.fac{Z.modes{pp}(1)};
                f_2 = sum(V(:));
                W = last_had{pp}.*G_transp_G{Z.modes{pp}(1)};
                f_3 = sum(W(:));
                fp(pp) = f_1 - 2* f_2 + f_3;
            else
                for kk=1:length(Z.size{Z.modes{pp}(2)})
                    fp(pp) = fp(pp) + norm(Z.object{pp}{kk}-G.fac{Z.modes{pp}(1)}*diag(G.fac{Z.modes{pp}(3)}(kk,:))*G.fac{Z.modes{pp}(2)}{kk}','fro')^2;
                end
            end
            fp(pp) = Z.weights(pp) *fp(pp);
        end
    end 
    fp
    f_tensors = sum(fp);
    
    if isfield(Z,'reg_func')
        for n = 1:nb_modes 
            if ~isempty(Z.reg_func{n})
                if iscell(G.constraint_fac{n}) %PAR2 2nd mode
                    if strcmp(Z.constraints{n}{1},'tPARAFAC2')
                        f_tensors = f_tensors + feval(Z.reg_func{n},G.fac{n});
                    else
                        for kk=1:length(G.constraint_fac{n})
                            f_tensors = f_tensors + feval(Z.reg_func{n},G.fac{n}{kk});
                        end
                    end
                else
                    f_tensors = f_tensors + feval(Z.reg_func{n},G.fac{n});
                end
            end
        end
    end
    
    if isfield(Z,'ridge')
         for n = 1:nb_modes
             if iscell(G.fac{n}) %PAR2 2nd mode
                    for kk=1:length(G.constraint_fac{n})
                        f_tensors = f_tensors + Z.ridge(n)*norm(G.fac{n}{kk},'fro')^2;
                    end
             else
                f_tensors = f_tensors + Z.ridge(n)*norm(G.fac{n},'fro')^2;
             end
         end
    end

    % residuals for coupling
    nb_couplings = max(Z.coupling.lin_coupled_modes);
    coupling_p = zeros(nb_couplings,1);
    for n = 1:nb_couplings
        ctype_n = Z.coupling.coupling_type(n);
        cmodes = find(Z.coupling.lin_coupled_modes==n);
        for jj = 1:length(cmodes)
            switch ctype_n
                case 0 %exact coupling
                    coupling_p(n) = coupling_p(n) + norm(G.fac{cmodes(jj)} - G.coupling_fac{n},'fro')/norm(G.fac{cmodes(jj)},'fro'); %no rho!
                case 1
                    coupling_p(n) = coupling_p(n) + norm(Z.coupling.coupl_trafo_matrices{cmodes(jj)}*G.fac{cmodes(jj)} - G.coupling_fac{n},'fro')/norm(Z.coupling.coupl_trafo_matrices{cmodes(jj)}*G.fac{cmodes(jj)},'fro'); %no rho!
                case 2
                    coupling_p(n) = coupling_p(n) + norm(G.fac{cmodes(jj)}*Z.coupling.coupl_trafo_matrices{cmodes(jj)} - G.coupling_fac{n},'fro')/norm(G.fac{cmodes(jj)}*Z.coupling.coupl_trafo_matrices{cmodes(jj)},'fro'); %no rho!
                case 3
                    coupling_p(n) = coupling_p(n) + norm(G.fac{cmodes(jj)} - Z.coupling.coupl_trafo_matrices{cmodes(jj)}*G.coupling_fac{n},'fro')/norm(G.fac{cmodes(jj)},'fro'); %no rho!
                case 4
                    coupling_p(n) = coupling_p(n) + norm(G.fac{cmodes(jj)} - G.coupling_fac{n}*Z.coupling.coupl_trafo_matrices{cmodes(jj)},'fro')/norm(G.fac{cmodes(jj)},'fro'); %no rho!
            end
        end
    end
    f_couplings = sum(coupling_p);
    if f_couplings>0
        f_couplings = f_couplings/length(find(coupling_p));
    end

    %residuals for constraints
    f_constraint_p = zeros(nb_modes,1);
    for n = 1:nb_modes 
        if ~isempty(G.constraint_fac{n})
            if iscell(G.constraint_fac{n}) %PAR2 2nd mode
                for kk=1:length(G.constraint_fac{n})
                    f_constraint_p(n) = f_constraint_p(n) + norm(G.fac{n}{kk}-G.constraint_fac{n}{kk},'fro')/norm(G.fac{n}{kk},'fro');
                end
                f_constraint_p(n) = f_constraint_p(n)/length(G.constraint_fac{n});
            else
                f_constraint_p(n) = norm(G.fac{n} - G.constraint_fac{n},'fro')/norm(G.fac{n},'fro');%no rho!
            end
        end
    end
    f_constraints = sum(f_constraint_p);
    if f_constraints>0
        f_constraints = f_constraints/length(find(f_constraint_p));
    end
    
    %residuals for Parafac2 internal coupling(s)
    f_PAR2_couplings_p = zeros(P,1);
    for pp = 1:P
        if strcmp(Z.model{pp},'PAR2')
            for kk=1:length(Z.size{Z.modes{pp}(2)})
                f_PAR2_couplings_p(pp) = f_PAR2_couplings_p(pp) + norm(G.fac{Z.modes{pp}(2)}{kk}-G.P{pp}{kk}*G.DeltaB{pp},'fro')/norm(G.fac{Z.modes{pp}(2)}{kk},'fro');
            end
        end
    end
    f_PAR2_couplings = sum(f_PAR2_couplings_p);
    if f_PAR2_couplings>0
        f_PAR2_couplings = f_PAR2_couplings/length(Z.size{Z.modes{pp}(2)}); %1/K
    end
end