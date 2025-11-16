function [inner_iter,lbfgsb_iterations,G] = ADMM_coupled_case1(Z,G,nb_modes,which_p,m,A,B,B2,coupled_modes,coupl_id,rho,options)
    inner_iter = 1;
    rel_primal_res_coupling = inf;
    rel_primal_res_constr = inf;
    rel_dual_res_coupling = inf;
    rel_dual_res_constr = inf;
    oldZ = cell(nb_modes,1);
    while (inner_iter<=options.MaxInnerIters &&(rel_primal_res_coupling>options.innerRelPrTol_coupl||rel_primal_res_constr>options.innerRelPrTol_constr||rel_dual_res_coupling>options.innerRelDualTol_coupl||rel_dual_res_constr>options.innerRelDualTol_constr))
        %exact coupling
        for mm=coupled_modes %update all factor matrices (can be done in parallel!)
            pp = which_p(mm);
            if strcmp(Z.loss_function{pp},'Frobenius')
                A_inner = A{mm} + rho{mm}/2*Z.coupling.coupl_trafo_matrices{mm}'*( G.coupling_fac{Z.coupling.lin_coupled_modes(mm)} - G.coupling_dual_fac{mm});
                if Z.constrained_modes(mm) %in case the mode is also constrained
                    A_inner = A_inner + rho{mm}/2*(G.constraint_fac{mm} - G.constraint_dual_fac{mm});
                end
                G.fac{mm} = sylvester(B2{mm},B{mm},A_inner); % solve Sylvester equation
                lbfgsb_iterations{m} = [];
            else
                [lbfgsb_iters(inner_iter)] = lbfgsb_update(pp,mm,Z.constrained_modes(mm),1,rho{mm}); %updates G.fac{m} with lbfgsb
                lbfgsb_iterations{mm} = lbfgsb_iters;
            end
        end
        
        % Update coupling factor (Delta) 
        oldDelta = G.coupling_fac{coupl_id};
        G.coupling_fac{coupl_id} = zeros(size(G.coupling_fac{coupl_id})); 
        sum_rho = 0;
        for jj = coupled_modes
            G.coupling_fac{coupl_id} = G.coupling_fac{coupl_id} + rho{jj}*(Z.coupling.coupl_trafo_matrices{jj}*G.fac{jj} + G.coupling_dual_fac{jj});
            sum_rho = sum_rho + rho{jj};
        end
        G.coupling_fac{coupl_id} = 1/sum_rho.*G.coupling_fac{coupl_id};
        
        % Update constraint factor (Z) and its dual (mu_Z) and mu_Delta
        for mm=coupled_modes % (can be done in parallel!)
            G.coupling_dual_fac{mm} = G.coupling_dual_fac{mm} + Z.coupling.coupl_trafo_matrices{mm}*G.fac{mm} - G.coupling_fac{Z.coupling.lin_coupled_modes(mm)}; % Update (mu_Delta)
            if Z.constrained_modes(mm)  
                oldZ{mm} = update_constraint(mm,rho{mm}); %updates G.constraint_fac{mm} and G.constraint_dual_fac{mm}
            end
        end
        inner_iter = inner_iter + 1; 
        [rel_primal_res_coupling,rel_dual_res_coupling] = eval_res_ADMM_coupl_case1(Z,G,coupled_modes,coupl_id,oldDelta);
        constrained_modes = coupled_modes(logical(Z.constrained_modes(coupled_modes)));
        if ~isempty(constrained_modes)
            [rel_primal_res_constr,rel_dual_res_constr] = eval_res_ADMM_constr(constrained_modes,oldZ); % does this work? integrate in loop instead? (no nested function)
        else
           rel_primal_res_constr = 0;
           rel_dual_res_constr =0;
        end
    end
    inner_iter = inner_iter-1;
end  