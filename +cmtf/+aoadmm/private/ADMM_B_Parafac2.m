function [inner_iter,G] = ADMM_B_Parafac2(Z,G,iter,A,L,m,p,rho,options)
%ADMM loop for mode B in Parafac2
% changes the global variable G 

    inner_iter = 1;
    rel_primal_res_constr = inf;
    rel_dual_res_constr = inf;
    rel_primal_res_coupling = inf;
    rel_dual_res_coupling = inf;
    oldP = cell(length(Z.size{Z.modes{p}(2)}),1);
    % ADMM loop
    while (inner_iter<=options.MaxInnerIters &&(rel_primal_res_coupling>options.innerRelPrTol_coupl||rel_primal_res_constr>options.innerRelPrTol_constr||rel_dual_res_coupling>options.innerRelDualTol_coupl||rel_dual_res_constr>options.innerRelDualTol_constr))
        rel_primal_res_constr = 0;
        rel_dual_res_constr = 0;
        rel_primal_res_coupling = 0;
        rel_dual_res_coupling = 0;
        for kk=1:length(Z.size{Z.modes{p}(2)})
            A_inner{kk} = A{kk} + rho(kk)/2*(G.P{p}{kk}*G.DeltaB{p} - G.mu_DeltaB{p}{kk});
            if Z.constrained_modes(m) && iter >= options.iter_start_PAR2Bkconstraint
                A_inner{kk} = A_inner{kk} + rho(kk)/2*(G.constraint_fac{m}{kk} - G.constraint_dual_fac{m}{kk});
            end
            G.fac{m}{kk} = (A_inner{kk}/L{kk}')/L{kk}; % forward-backward substitution
            % update P_k
            [U,~,V] = svd((G.fac{m}{kk} + G.mu_DeltaB{p}{kk})*G.DeltaB{p}','econ');
            oldP{kk} = G.P{p}{kk};
            G.P{p}{kk} = U*V';
        end
        % update DeltaB
        oldDeltaB = G.DeltaB{p};
        G.DeltaB{p} = zeros(size(G.DeltaB{p}));
        sum_rho_k = 0;
        for kk=1:length(Z.size{Z.modes{p}(2)})
            G.DeltaB{p} = G.DeltaB{p} + rho(kk)*G.P{p}{kk}'*(G.fac{m}{kk} + G.mu_DeltaB{p}{kk});
            sum_rho_k = sum_rho_k + rho(kk);
        end
        G.DeltaB{p} = G.DeltaB{p}./sum_rho_k;
        for kk=1:length(Z.size{Z.modes{p}(2)})
            G.mu_DeltaB{p}{kk} = G.mu_DeltaB{p}{kk} + G.fac{m}{kk} - G.P{p}{kk}*G.DeltaB{p};
        end
        
        % Update constraint factor (Z_B) and its dual (mu_Z_B) if
        % constrained
        if Z.constrained_modes(m) && iter >= options.iter_start_PAR2Bkconstraint
            oldZ = G.constraint_fac{m};
            if strcmp(Z.constraints{m}{1},'tPARAFAC2')
                G.constraint_fac{m} = feval(Z.prox_operators{m},cellfun(@(x, y) x + y, G.fac{m}, G.constraint_dual_fac{m}, 'UniformOutput', false),rho);
                for kk=1:length(Z.size{Z.modes{p}(2)})
                    G.constraint_dual_fac{m}{kk} = G.constraint_dual_fac{m}{kk} + G.fac{m}{kk} - G.constraint_fac{m}{kk};             
                    % sum up residuals
                    rel_primal_res_constr = rel_primal_res_constr + norm(G.fac{m}{kk} - G.constraint_fac{m}{kk},'fro')/norm(G.fac{m}{kk},'fro')/length(Z.size{Z.modes{p}(2)});
                    scaling = norm(G.constraint_dual_fac{m}{kk},'fro');
                    if scaling>0
                        rel_dual_res_constr = rel_dual_res_constr + norm(oldZ{kk} - G.constraint_fac{m}{kk},'fro')/scaling/length(Z.size{Z.modes{p}(2)});
                    else
                        rel_dual_res_constr = rel_dual_res_constr + norm(oldZ{kk} - G.constraint_fac{m}{kk},'fro')/length(Z.size{Z.modes{p}(2)});
                    end
                end
            else
                for kk=1:length(Z.size{Z.modes{p}(2)})
                    G.constraint_fac{m}{kk} = feval(Z.prox_operators{m},(G.fac{m}{kk} + G.constraint_dual_fac{m}{kk}),rho(kk));
                    G.constraint_dual_fac{m}{kk} = G.constraint_dual_fac{m}{kk} + G.fac{m}{kk} - G.constraint_fac{m}{kk};             
                    % sum up residuals
                    rel_primal_res_constr = rel_primal_res_constr + norm(G.fac{m}{kk} - G.constraint_fac{m}{kk},'fro')/norm(G.fac{m}{kk},'fro')/length(Z.size{Z.modes{p}(2)});
                    scaling = norm(G.constraint_dual_fac{m}{kk},'fro');
                    if scaling>0
                        rel_dual_res_constr = rel_dual_res_constr + norm(oldZ{kk} - G.constraint_fac{m}{kk},'fro')/scaling/length(Z.size{Z.modes{p}(2)});
                    else
                        rel_dual_res_constr = rel_dual_res_constr + norm(oldZ{kk} - G.constraint_fac{m}{kk},'fro')/length(Z.size{Z.modes{p}(2)});
                    end
                end
            end
        end
        
        for kk=1:length(Z.size{Z.modes{p}(2)})
            rel_primal_res_coupling = rel_primal_res_coupling + norm(G.fac{m}{kk} - G.P{p}{kk}*G.DeltaB{p},'fro')/norm(G.fac{m}{kk},'fro')/length(Z.size{Z.modes{p}(2)});
            rel_dual_res_coupling = rel_dual_res_coupling + norm(oldP{kk}*oldDeltaB - G.P{p}{kk}*G.DeltaB{p},'fro')/norm(G.mu_DeltaB{p}{kk},'fro')/length(Z.size{Z.modes{p}(2)});
        end
        inner_iter = inner_iter + 1; 
    end
    inner_iter = inner_iter-1;
end