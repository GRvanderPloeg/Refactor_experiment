function [inner_iter,lbfgsb_iterations,G] = ADMM_coupled(Z,G,nb_modes,which_p,ctype,lscalar,uscalar,fh,gh,A,B,coupled_modes,coupl_id,rho,options)
%ADMM loop for linearly coupled modes, all coupling types (0-4).
%
% Coupling types:
%   0: C = Delta          (exact coupling)
%   1: H*C = Delta
%   2: C*H = Delta
%   3: C = H*Delta
%   4: C = Delta*H
%
% B is the pre-coupling Gram-matrix term (w*C + ridge + bsum) for each
% coupled mode.  ADMM_coupled adds the coupling-specific rho terms and
% computes the Cholesky factor L (or B2 for the Sylvester case) internally,
% once before the inner loop.

    inner_iter        = 1;
    rel_primal_res_coupling = inf;
    rel_primal_res_constr   = inf;
    rel_dual_res_coupling   = inf;
    rel_dual_res_constr     = inf;
    oldZ              = cell(nb_modes,1);
    lbfgsb_iterations = cell(nb_modes,1);

    %% Precompute L (Cholesky) or B2 (Sylvester) — done once before the loop
    L  = cell(nb_modes,1);
    B2 = cell(nb_modes,1);
    switch ctype
        case 0  % exact coupling: B += rho/2*I (twice if also constrained)
            for mm = coupled_modes
                pp = which_p(mm);
                if strcmp(Z.loss_function{pp},'Frobenius')
                    if strcmp(Z.model{pp},'PAR2') && 3 == find(Z.modes{pp}==mm)
                        for kk = 1:length(Z.size{Z.modes{pp}(2)})
                            B{mm}{kk} = B{mm}{kk} + rho{mm}(kk)/2*eye(size(B{mm}{kk}));
                            if Z.constrained_modes(mm)
                                B{mm}{kk} = B{mm}{kk} + rho{mm}(kk)/2*eye(size(B{mm}{kk}));
                            end
                            L{mm}{kk} = chol(B{mm}{kk}','lower');
                        end
                    else
                        B{mm} = B{mm} + rho{mm}/2*eye(size(B{mm}));
                        if Z.constrained_modes(mm)
                            B{mm} = B{mm} + rho{mm}/2*eye(size(B{mm}));
                        end
                        L{mm} = chol(B{mm}','lower');
                    end
                end
            end
        case 1  % HC=Delta: Sylvester equation; B2 = rho/2*H'*H [+ constraint term]
            for mm = coupled_modes
                pp = which_p(mm);
                if strcmp(Z.loss_function{pp},'Frobenius')
                    B2{mm} = rho{mm}/2*Z.coupling.coupl_trafo_matrices{mm}'*Z.coupling.coupl_trafo_matrices{mm};
                    if Z.constrained_modes(mm)
                        B2{mm} = B2{mm} + rho{mm}/2*eye(size(B2{mm}));
                    end
                end
            end
        case 2  % CH=Delta: B += rho/2*H*H'
            for mm = coupled_modes
                pp = which_p(mm);
                if strcmp(Z.loss_function{pp},'Frobenius')
                    B{mm} = B{mm} + rho{mm}/2*Z.coupling.coupl_trafo_matrices{mm}*Z.coupling.coupl_trafo_matrices{mm}';
                    if Z.constrained_modes(mm)
                        B{mm} = B{mm} + rho{mm}/2*eye(size(B{mm}));
                    end
                    L{mm} = chol(B{mm}','lower');
                end
            end
        case {3,4}  % C=HDelta or C=DeltaH: B += rho/2*I
            for mm = coupled_modes
                pp = which_p(mm);
                if strcmp(Z.loss_function{pp},'Frobenius')
                    B{mm} = B{mm} + rho{mm}/2*eye(size(B{mm}));
                    if Z.constrained_modes(mm)
                        B{mm} = B{mm} + rho{mm}/2*eye(size(B{mm}));
                    end
                    L{mm} = chol(B{mm}','lower');
                end
            end
    end

    %% Inner ADMM loop
    while (inner_iter <= options.MaxInnerIters && ...
          (rel_primal_res_coupling > options.innerRelPrTol_coupl  || ...
           rel_primal_res_constr   > options.innerRelPrTol_constr || ...
           rel_dual_res_coupling   > options.innerRelDualTol_coupl || ...
           rel_dual_res_constr     > options.innerRelDualTol_constr))

        %% Step 1: Update factor matrices
        for mm = coupled_modes
            pp = which_p(mm);
            if strcmp(Z.loss_function{pp},'Frobenius')
                n = Z.coupling.lin_coupled_modes(mm);
                switch ctype
                    case 0  % exact coupling
                        if strcmp(Z.model{pp},'PAR2') && 3 == find(Z.modes{pp}==mm)
                            for kk = 1:length(Z.size{Z.modes{pp}(2)})
                                A_inner = A{mm}{kk} + rho{mm}(kk)/2*(G.coupling_fac{n}(kk,:)' - G.coupling_dual_fac{mm}(kk,:)');
                                if Z.constrained_modes(mm)
                                    A_inner = A_inner + rho{mm}(kk)/2*(G.constraint_fac{mm}(kk,:)' - G.constraint_dual_fac{mm}(kk,:)');
                                end
                                G.fac{mm}(kk,:) = (L{mm}{kk}'\(L{mm}{kk}\A_inner))';
                            end
                        else
                            A_inner = A{mm} + rho{mm}/2*(G.coupling_fac{n} - G.coupling_dual_fac{mm});
                            if Z.constrained_modes(mm)
                                A_inner = A_inner + rho{mm}/2*(G.constraint_fac{mm} - G.constraint_dual_fac{mm});
                            end
                            G.fac{mm} = (A_inner/L{mm}')/L{mm};
                        end
                    case 1  % HC=Delta: Sylvester
                        A_inner = A{mm} + rho{mm}/2*Z.coupling.coupl_trafo_matrices{mm}'*(G.coupling_fac{n} - G.coupling_dual_fac{mm});
                        if Z.constrained_modes(mm)
                            A_inner = A_inner + rho{mm}/2*(G.constraint_fac{mm} - G.constraint_dual_fac{mm});
                        end
                        G.fac{mm} = sylvester(B2{mm},B{mm},A_inner);
                    case 2  % CH=Delta
                        A_inner = A{mm} + rho{mm}/2*(G.coupling_fac{n} - G.coupling_dual_fac{mm})*Z.coupling.coupl_trafo_matrices{mm}';
                        if Z.constrained_modes(mm)
                            A_inner = A_inner + rho{mm}/2*(G.constraint_fac{mm} - G.constraint_dual_fac{mm});
                        end
                        G.fac{mm} = (A_inner/L{mm}')/L{mm};
                    case 3  % C=HDelta
                        A_inner = A{mm} + rho{mm}/2*(Z.coupling.coupl_trafo_matrices{mm}*G.coupling_fac{n} - G.coupling_dual_fac{mm});
                        if Z.constrained_modes(mm)
                            A_inner = A_inner + rho{mm}/2*(G.constraint_fac{mm} - G.constraint_dual_fac{mm});
                        end
                        G.fac{mm} = (A_inner/L{mm}')/L{mm};
                    case 4  % C=DeltaH
                        A_inner = A{mm} + rho{mm}/2*(G.coupling_fac{n}*Z.coupling.coupl_trafo_matrices{mm} - G.coupling_dual_fac{mm});
                        if Z.constrained_modes(mm)
                            A_inner = A_inner + rho{mm}/2*(G.constraint_fac{mm} - G.constraint_dual_fac{mm});
                        end
                        G.fac{mm} = (A_inner/L{mm}')/L{mm};
                end
                lbfgsb_iterations{mm} = [];
            else  % non-Frobenius: L-BFGS-B
                [lbfgsb_iters(inner_iter),G] = lbfgsb_update(Z,G,lscalar,uscalar,fh,gh,pp,mm,Z.constrained_modes(mm),ctype,rho{mm},options);
                lbfgsb_iterations{mm} = lbfgsb_iters;
            end
        end

        %% Step 2: Update coupling variable (Delta)
        oldDelta = G.coupling_fac{coupl_id};
        switch ctype
            case 0  % Delta = mean_rho(C + mu)  [PAR2 third-mode variant included]
                G.coupling_fac{coupl_id} = zeros(size(G.coupling_fac{coupl_id}));
                sum_rho = 0;
                for jj = coupled_modes
                    pp = which_p(jj);
                    if strcmp(Z.model{pp},'PAR2') && 3 == find(Z.modes{pp}==jj)
                        for kk = 1:length(Z.size{Z.modes{pp}(2)})
                            G.coupling_fac{coupl_id}(kk,:) = G.coupling_fac{coupl_id}(kk,:) + rho{jj}(kk)*(G.fac{jj}(kk,:) + G.coupling_dual_fac{jj}(kk,:));
                        end
                    else
                        G.coupling_fac{coupl_id} = G.coupling_fac{coupl_id} + rho{jj}*(G.fac{jj} + G.coupling_dual_fac{jj});
                    end
                    sum_rho = sum_rho + rho{jj};
                end
                G.coupling_fac{coupl_id} = 1./sum_rho'.*G.coupling_fac{coupl_id};
            case {1,2}  % Delta = mean_rho(H*C + mu) or mean_rho(C*H + mu)
                G.coupling_fac{coupl_id} = zeros(size(G.coupling_fac{coupl_id}));
                sum_rho = 0;
                for jj = coupled_modes
                    if ctype == 1
                        contrib = Z.coupling.coupl_trafo_matrices{jj}*G.fac{jj} + G.coupling_dual_fac{jj};
                    else
                        contrib = G.fac{jj}*Z.coupling.coupl_trafo_matrices{jj} + G.coupling_dual_fac{jj};
                    end
                    G.coupling_fac{coupl_id} = G.coupling_fac{coupl_id} + rho{jj}*contrib;
                    sum_rho = sum_rho + rho{jj};
                end
                G.coupling_fac{coupl_id} = 1/sum_rho.*G.coupling_fac{coupl_id};
            case 3  % Delta = (sum rho*H'*H) \ (sum rho*H'*(C + mu))
                AA = zeros(size(Z.coupling.coupl_trafo_matrices{coupled_modes(1)},2));
                BB = zeros(size(Z.coupling.coupl_trafo_matrices{coupled_modes(1)},2), size(G.fac{coupled_modes(1)},2));
                for jj = coupled_modes
                    AA = AA + rho{jj}*Z.coupling.coupl_trafo_matrices{jj}'*Z.coupling.coupl_trafo_matrices{jj};
                    BB = BB + rho{jj}*Z.coupling.coupl_trafo_matrices{jj}'*(G.fac{jj} + G.coupling_dual_fac{jj});
                end
                G.coupling_fac{coupl_id} = AA\BB;
            case 4  % Delta = (sum rho*(C + mu)*H') / (sum rho*H*H')
                AA = zeros(size(Z.coupling.coupl_trafo_matrices{coupled_modes(1)},1));
                BB = zeros(size(G.fac{coupled_modes(1)},1), size(Z.coupling.coupl_trafo_matrices{coupled_modes(1)},1));
                for jj = coupled_modes
                    AA = AA + rho{jj}*Z.coupling.coupl_trafo_matrices{jj}*Z.coupling.coupl_trafo_matrices{jj}';
                    BB = BB + rho{jj}*(G.fac{jj} + G.coupling_dual_fac{jj})*Z.coupling.coupl_trafo_matrices{jj}';
                end
                G.coupling_fac{coupl_id} = BB/AA;
        end

        %% Step 3: Update coupling dual (mu_Delta) and constraint variables
        for mm = coupled_modes
            n = Z.coupling.lin_coupled_modes(mm);
            switch ctype
                case 0
                    G.coupling_dual_fac{mm} = G.coupling_dual_fac{mm} + G.fac{mm} - G.coupling_fac{n};
                case 1
                    G.coupling_dual_fac{mm} = G.coupling_dual_fac{mm} + Z.coupling.coupl_trafo_matrices{mm}*G.fac{mm} - G.coupling_fac{n};
                case 2
                    G.coupling_dual_fac{mm} = G.coupling_dual_fac{mm} + G.fac{mm}*Z.coupling.coupl_trafo_matrices{mm} - G.coupling_fac{n};
                case 3
                    G.coupling_dual_fac{mm} = G.coupling_dual_fac{mm} + G.fac{mm} - Z.coupling.coupl_trafo_matrices{mm}*G.coupling_fac{n};
                case 4
                    G.coupling_dual_fac{mm} = G.coupling_dual_fac{mm} + G.fac{mm} - G.coupling_fac{n}*Z.coupling.coupl_trafo_matrices{mm};
            end
            if Z.constrained_modes(mm)
                [oldZ{mm},G] = update_constraint(Z,G,mm,rho{mm});
            end
        end

        %% Step 4: Evaluate convergence
        inner_iter = inner_iter + 1;
        [rel_primal_res_coupling,rel_dual_res_coupling] = eval_res_ADMM_coupl(Z,G,coupled_modes,coupl_id,ctype,oldDelta);
        constrained_modes = coupled_modes(logical(Z.constrained_modes(coupled_modes)));
        if ~isempty(constrained_modes)
            [rel_primal_res_constr,rel_dual_res_constr] = eval_res_ADMM_constr(G,constrained_modes,oldZ);
        else
            rel_primal_res_constr = 0;
            rel_dual_res_constr   = 0;
        end
    end
    inner_iter = inner_iter - 1;
end
