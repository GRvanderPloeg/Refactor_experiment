function [inner_iter,lbfgsb_iterations,G] = ADMM_coupled_case0(Z,G,nb_modes,which_p,m,lscalar,uscalar,fh,gh,A,L,coupled_modes,coupl_id,rho,options)
    inner_iter = 1;
    rel_primal_res_coupling = inf;
    rel_primal_res_constr = inf;
    rel_dual_res_coupling = inf;
    rel_dual_res_constr = inf;
    oldZ = cell(nb_modes,1);
    lbfgsb_iterations = cell(nb_modes,1);
    while (inner_iter<=options.MaxInnerIters &&(rel_primal_res_coupling>options.innerRelPrTol_coupl||rel_primal_res_constr>options.innerRelPrTol_constr||rel_dual_res_coupling>options.innerRelDualTol_coupl||rel_dual_res_constr>options.innerRelDualTol_constr))
        %exact coupling
        for mm=coupled_modes %update all factor matrices (can be done in parallel!)
            pp = which_p(mm);
            if strcmp(Z.loss_function{pp},'Frobenius')
                if strcmp(Z.model{pp},'PAR2')  && 3 == find(Z.modes{pp}==mm) %third Parafac2 mode
                    for kk=1:length(Z.size{Z.modes{pp}(2)})
                        A_inner = A{mm}{kk} + rho{mm}(kk)/2*( G.coupling_fac{Z.coupling.lin_coupled_modes(mm)}(kk,:)' - G.coupling_dual_fac{mm}(kk,:)');
                        if Z.constrained_modes(mm) %in case the mode is also constrained
                             A_inner = A_inner + rho{mm}(kk)/2*(G.constraint_fac{mm}(kk,:)' - G.constraint_dual_fac{mm}(kk,:)');
                        end
                    G.fac{mm}(kk,:) = (L{mm}{kk}'\(L{mm}{kk}\A_inner))'; % forward-backward substitution
                    end
                else
                    A_inner = A{mm} + rho{mm}/2*( G.coupling_fac{Z.coupling.lin_coupled_modes(mm)} - G.coupling_dual_fac{mm});
                    if Z.constrained_modes(mm) %in case the mode is also constrained
                         A_inner = A_inner+ rho{mm}/2*(G.constraint_fac{mm} - G.constraint_dual_fac{mm});
                    end
                    G.fac{mm} = (A_inner/L{mm}')/L{mm}; % forward-backward substitution
                end
                lbfgsb_iterations{m} = [];
            else
                [lbfgsb_iters(inner_iter),G] = lbfgsb_update(Z,G,lscalar,uscalar,fh,gh,pp,mm,Z.constrained_modes(mm),0,rho{mm},options); %updates G.fac{m}
                lbfgsb_iterations{mm} = lbfgsb_iters;
            end
        end
        
        % Update coupling factor (Delta) 
        oldDelta = G.coupling_fac{coupl_id};
        G.coupling_fac{coupl_id} = zeros(size(G.coupling_fac{coupl_id})); 
        sum_rho = 0;
        for jj = coupled_modes
            pp = which_p(jj);
            if strcmp(Z.model{pp},'PAR2')  && 3 == find(Z.modes{pp}==jj) %third Parafac2 mode
                for kk=1:length(Z.size{Z.modes{pp}(2)})
                    G.coupling_fac{coupl_id}(kk,:) = G.coupling_fac{coupl_id}(kk,:) + rho{jj}(kk)*(G.fac{jj}(kk,:) + G.coupling_dual_fac{jj}(kk,:));
                end
            else
                G.coupling_fac{coupl_id} = G.coupling_fac{coupl_id} + rho{jj}*(G.fac{jj} + G.coupling_dual_fac{jj});
            end
            sum_rho = sum_rho + rho{jj};
        end
        G.coupling_fac{coupl_id} = 1./sum_rho'.*G.coupling_fac{coupl_id};
        
        % Update constraint factor (Z) and its dual (mu_Z) and mu_Delta
        for mm=coupled_modes % (can be done in parallel!)
            G.coupling_dual_fac{mm} = G.coupling_dual_fac{mm} + G.fac{mm} - G.coupling_fac{Z.coupling.lin_coupled_modes(mm)}; % Update (mu_Delta)
            if Z.constrained_modes(mm)  
                [oldZ{mm}, G] = update_constraint(Z,G,mm,rho{mm}); %updates G.constraint_fac{mm} and G.constraint_dual_fac{mm}
            end
        end
        inner_iter = inner_iter + 1; 
        [rel_primal_res_coupling,rel_dual_res_coupling] = eval_res_ADMM_coupl_case0(G,coupled_modes,coupl_id,oldDelta);
        constrained_modes = coupled_modes(logical(Z.constrained_modes(coupled_modes)));
        if ~isempty(constrained_modes)
            [rel_primal_res_constr,rel_dual_res_constr] = eval_res_ADMM_constr(G,constrained_modes,oldZ); % does this work? integrate in loop instead? (no nested function)
        else
           rel_primal_res_constr = 0;
           rel_dual_res_constr =0;
        end
    end
    inner_iter = inner_iter-1;
end  