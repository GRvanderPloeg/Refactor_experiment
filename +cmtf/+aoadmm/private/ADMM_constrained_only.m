function [inner_iter,lbfgsb_iterations,G] = ADMM_constrained_only(Z,G,nb_modes,A,L,m,p,rho,options)
%ADMM loop for mode m, where mode m is constrained, but not coupled!
% changes the global variable G (only fields related to mode m)

    inner_iter = 1;
    rel_primal_res_constr = inf;
    rel_dual_res_constr = inf;
    oldZ = cell(nb_modes,1);
    % ADMM loop
    while (inner_iter<=options.MaxInnerIters &&(rel_primal_res_constr>options.innerRelPrTol_constr||rel_dual_res_constr>options.innerRelDualTol_constr))
        if strcmp(Z.loss_function{p},'Frobenius')
            if strcmp(Z.model{p},'PAR2')  && 3 == find(Z.modes{p}==m) %third Parafac2 mode
                for kk=1:length(Z.size{Z.modes{p}(2)})
                    A_inner = A{kk} + rho{m}(kk)/2*(G.constraint_fac{m}(kk,:)' - G.constraint_dual_fac{m}(kk,:)');
                    G.fac{m}(kk,:) = (L{kk}'\(L{kk}\A_inner))'; % forward-backward substitution
                end
            else
                A_inner = A + rho{m}/2*(G.constraint_fac{m} - G.constraint_dual_fac{m});
                G.fac{m} = (A_inner/L')/L; % forward-backward substitution
            end
            lbfgsb_iterations{m} = [];
        else % other loss function, use lbfgsb
            [lbfgsb_iterations{m}(inner_iter)] = lbfgsb_update(p,m,true,-1,rho{m}); %updates G.fac{m}
        end

        % Update constraint factor (Z) and its dual (mu_Z)
        [oldZ{m},G] = update_constraint(Z,G,m,rho{m}); %updates G.constraint_fac{mm} and G.constraint_dual_fac{mm}

        inner_iter = inner_iter + 1; 
        [rel_primal_res_constr,rel_dual_res_constr] = eval_res_ADMM_constr(G,m,oldZ); 
    end
    inner_iter = inner_iter-1;
end