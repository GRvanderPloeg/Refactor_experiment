function [G,out] = cmtf_fun_AOADMM(Z,Znorm_const, G,fh,gh,lscalar,uscalar,options)
% Alternating Optimization ADMM (AOADMM) inner loop for CMTF.
%
% Performs one full pass of alternating optimization over all modes of the
% coupled tensor factorization problem, using ADMM to handle constraints
% and couplings. Supports CP and PARAFAC2 models, multiple coupling types,
% Frobenius and non-Frobenius loss functions.
%
% Syntax:
%   [G, out] = cmtf_fun_AOADMM(Z, Znorm_const, G, fh, gh, lscalar, uscalar, options)
%
% Inputs:
%   Z          - Problem structure with fields: object, modes, size, weights,
%                coupling, loss_function, model, constrained_modes, (optional) ridge
%   Znorm_const - Scalar normalization constant used in objective evaluation
%   G          - Factor struct with field .fac (cell array of factor matrices)
%   fh         - Function handle for non-Frobenius loss evaluation
%   gh         - Gradient handle for non-Frobenius loss (used by L-BFGS-B)
%   lscalar    - Lower bound scalar for L-BFGS-B updates
%   uscalar    - Upper bound scalar for L-BFGS-B updates
%   options    - Struct with solver settings (MaxOuterIters, Display,
%                DisplayIters, AbsFuncTol, OuterRelTol, bsum, etc.)
%
% Outputs:
%   G   - Updated factor struct after all AO-ADMM iterations
%   out - Struct with convergence diagnostics: f_tensors, f_couplings,
%         f_constraints, f_PAR2_couplings, exit_flag, OuterIterations,
%         func_val_conv, innerIters, time_at_it, (optionally lbfgsb_iterations)

    if isfield(options,'lbfgsb_options')
        lbfgsb_options = options.lbfgsb_options;
    end
    if ~isfield(options,'iter_start_PAR2Bkconstraint')
        options.iter_start_PAR2Bkconstraint = 0;
    end
    couplings = unique(Z.coupling.lin_coupled_modes);
    nb_modes  = length(Z.size);
    which_p = zeros(1,nb_modes);
    for i=1:nb_modes
        which_p(i) = find(cellfun(@(x) any(ismember(x,i)),Z.modes));
    end
    P = numel(Z.object);
    G_transp_G = cell(nb_modes,1);
    sum_column_norms_sqr = zeros(nb_modes,1);
    A = cell(nb_modes,1);
    C = cell(nb_modes,1);
    B = cell(nb_modes,1);
    B2 = cell(nb_modes,1);
    L = cell(nb_modes,1);
    rho = cell(nb_modes,1);
    last_m = zeros(P,1);
    last_mttkrp = cell(P,1); %for efficient function value computation
    last_had = cell(P,1); %for efficient function value computation
    out.innerIters = zeros(nb_modes,1);

    [f_tensors,f_couplings,f_constraints,f_PAR2_couplings] = CMTF_AOADMM_func_eval(Z,G,G_transp_G,Znorm_const,[],[],[],fh);
    f_total = f_tensors+f_couplings+f_constraints + f_PAR2_couplings;
    func_val(1) = f_tensors;
    func_coupl(1) = f_couplings;
    func_constr(1) = f_constraints;
    func_PAR2_coupl(1) = f_PAR2_couplings;
    tstart = tic;
    time_at_it(1) = 0;
    %display first iteration
    if strcmp(options.Display,'iter') || strcmp(options.Display,'final')
        fprintf(1,' Iter  f total      f tensors      f couplings    f constraints    f PAR2 couplings\n');
        fprintf(1,'------ ------------ -------------  -------------- ---------------- ----------------\n');
    end

    if strcmp(options.Display,'iter')
        fprintf(1,'%6d %12f %12f %12f %17f %12f\n', 0, f_total, f_tensors, f_couplings,f_constraints, f_PAR2_couplings);
    end
    iter = 1;
    
    for m=1:nb_modes
        p = which_p(m);
        if strcmp(Z.loss_function{p},'Frobenius')
            if strcmp(Z.model{p},'CP')
                G_transp_G{m} = G.fac{m}'*G.fac{m}; %precompute G'*G;
            elseif strcmp(Z.model{p},'PAR2')
                if 1 == find(Z.modes{p}==m)
                    G_transp_G{m} = G.fac{m}'*G.fac{m};
                elseif 2 == find(Z.modes{p}==m)
                    for k=1:length(Z.size{m})
                        G_transp_G{m}{k} = G.fac{m}{k}'*G.fac{m}{k};
                    end
                end
            end
        else
            for r=1:size(G.fac{m},2)
                sum_column_norms_sqr(m,1) = sum_column_norms_sqr(m,1)+norm(G.fac{m}(:,r))^2; % need this to compute rho in the non-Frobenius case
            end
        end
    end

    s = warning('error', 'MATLAB:nearlySingularMatrix');
    illconditioned = 0;
    
    stop = false;
    while(iter<=options.MaxOuterIters && ~stop)    
        %try
        for coupl_id=couplings %loop over all couplings (and non-coupled modes, if couplings=0)
            coupled_modes = find(Z.coupling.lin_coupled_modes==coupl_id); % all modes with this coupling_id
            for p=unique(which_p(coupled_modes)) % loop over all coupled tensors for this coupling_id (can be done in parallel!)
                if strcmp(Z.model{p},'CP')

                    % AO: mode 1 (the coupled mode)
                    for m=coupled_modes(which_p(coupled_modes)==p) %loop over all modes in tensor p with this coupling_id (can NOT be done in parallel)
                        if strcmp(Z.loss_function{p},'Frobenius')
                            % precomputations
                            if length(size(Z.object{p}))>=3  % Tensor
                                A{m} = Z.weights(p) *mttkrp(Z.object{p},G.fac(Z.modes{p}),find(Z.modes{p}==m)); %efficient calculation of matricized tensor with kathrirao product of all factor matrices, but the mth
                                C{m} = ones(size(G_transp_G{m}));
                                for j=Z.modes{p} 
                                    if(j~=m)
                                        C{m} = C{m} .* G_transp_G{j}; % efficient calculation of the product of khatrirao products
                                    end
                                end

                            else % Matrix
                                matrix_mode = find(Z.modes{p}==m);
                                if matrix_mode == 1 %first mode in matrix
                                    A{m} = Z.weights(p)* double(Z.object{p})*G.fac{Z.modes{p}(2)};
                                    C{m} = G_transp_G{Z.modes{p}(2)};
                                else %second mode in matrix
                                    A{m} = Z.weights(p)* double(Z.object{p})'*G.fac{Z.modes{p}(1)}; %transposed M!
                                    C{m} = G_transp_G{Z.modes{p}(1)};
                                end
                            end
                            rho{m} = trace(C{m})/size(C{m},1);
                            B{m} = Z.weights(p)* C{m}; 
                            if isfield(Z,'ridge')
                                B{m} = B{m} + Z.ridge(m)*eye(size(B{m}));
                            end

                            last_mttkrp{p} = A{m}*1/Z.weights(p);
                            last_had{p} = C{m};
                            last_m(p) = m;
                            if options.bsum
                                A{m} = A{m} + options.bsum_weight/2*G.fac{m};
                                B{m} = B{m} + options.bsum_weight/2*eye(size(B{m}));
                            end
                        else % other loss than Frobenius
                            rho{m} = sum(sum_column_norms_sqr([1:m-1,m+1:end])); 
                        end
                        if coupl_id==0 %modes are not coupled
                            if (Z.constrained_modes(m)==0) % mode is not constrained
                                if strcmp(Z.loss_function{p},'Frobenius')
                                    G.fac{m} = A{m}/B{m}; % Update factor matrices (Least squares update)
                                else
                                    [lbfgsb_iterations,G] = lbfgsb_update(Z,G,lscalar,uscalar,fh,gh,p,m,false,-1,rho{m},options); %updates G.fac{m} with lbfgsb
                                end
                                inner_iters = 1;
                            else % mode is constrained, use ADMM
                                if strcmp(Z.loss_function{p},'Frobenius')
                                    B{m} = B{m}+ rho{m}/2*eye(size(B{m})); % for constraint
                                    L{m} = chol(B{m}','lower'); %precompute Cholesky decomposition of B (only works in the chase when rho does not change between inner iterations)
                                end
                                [inner_iters,lbfgsb_iterations,G] = ADMM_constrained_only(Z,G,nb_modes,A{m},L{m},m,p,rho,options);
                            end
                            out.innerIters(m,iter)= inner_iters;
                            if strcmp(Z.loss_function{p},'Frobenius')
                                G_transp_G{m} = G.fac{m}'*G.fac{m}; % update G transposed G for mth mode
                            else
                                out.lbfgsb_iterations{m,iter} = lbfgsb_iterations;
                                for r=1:size(G.fac{m},2)
                                    sum_column_norms_sqr(m,1) = norm(G.fac{m}(:,r))^2;
                                end
                            end
                        end
                    end
                elseif strcmp(Z.model{p},'PAR2')
                    for m=coupled_modes(which_p(coupled_modes)==p)
                        if 1 == find(Z.modes{p}==m) %first Parafac2 mode (can be coupled)
                           A{m} = zeros(size(G.fac{m}));
                           C{m} = zeros(size(G.fac{m},2),size(G.fac{m},2));
                           for k=1:length(Z.size{Z.modes{p}(2)})
                               A{m} = A{m} + Z.object{p}{k}*G.fac{Z.modes{p}(2)}{k}*diag(G.fac{Z.modes{p}(3)}(k,:));
                               C{m} = C{m} + diag(G.fac{Z.modes{p}(3)}(k,:))*G_transp_G{Z.modes{p}(2)}{k}*diag(G.fac{Z.modes{p}(3)}(k,:));
                           end 
                           last_had{p} = C{m};
                           last_mttkrp{p} = A{m};
                           last_m(p) = 1;
                           A{m} = Z.weights(p)*A{m};
                           rho{m} = trace(C{m})/size(C{m},1);
                           B{m} = Z.weights(p)* C{m};
                           if isfield(Z,'ridge')
                                B{m} = B{m} + Z.ridge(m)*eye(size(B{m}));
                            end
                           if options.bsum
                                A{m} = A{m} + options.bsum_weight/2*G.fac{m};
                                B{m} = B{m} + options.bsum_weight/2*eye(size(B{m}));
                           end
                           if coupl_id==0 %mode is not coupled
                                if (Z.constrained_modes(m)==0) % mode is not constrained
                                    G.fac{m} = A{m}/B{m}; %ALS update
                                    inner_iters = 1;
                                else % mode is constrained, use ADMM
                                    B{m} = B{m}+ rho{m}/2*eye(size(B{m})); % for constraint
                                    L{m} = chol(B{m}','lower'); %precompute Cholesky decomposition
                                    [inner_iters,lbfgsb_iterations,G] = ADMM_constrained_only(Z,G,nb_modes,A{m},L{m},m,p,rho,options);
                                end
                            end
                            out.innerIters(m,iter)= inner_iters;
                            G_transp_G{m} = G.fac{m}'*G.fac{m}; % update G transposed G for mth mode
                        elseif 2 == find(Z.modes{p}==m) %second Parafac2 mode (the funny mode) no external coupling allowed
                            for k=1:length(Z.size{Z.modes{p}(2)})
                                A{m}{k} = Z.weights(p)*Z.object{p}{k}'*G.fac{Z.modes{p}(1)}*diag(G.fac{Z.modes{p}(3)}(k,:));
                                C{m}{k} = diag(G.fac{Z.modes{p}(3)}(k,:))*G_transp_G{Z.modes{p}(1)}*diag(G.fac{Z.modes{p}(3)}(k,:));
                                rho{m}(k) = trace(C{m}{k})/size(C{m}{k},1);
                                if isfield(options, 'increase_factor_rhoBk')
                                    rho{m}(k) = options.increase_factor_rhoBk * rho{m}(k);
                                end
                                B{m}{k} = Z.weights(p)* C{m}{k};
                                B{m}{k} = B{m}{k} + rho{m}(k)/2*eye(size(B{m}{k})); %always coupled 
                                if isfield(Z,'ridge')
                                    B{m}{k} = B{m}{k} + Z.ridge(m)*eye(size(B{m}{k}));
                                end
                                if options.bsum
                                    A{m}{k} = A{m}{k} + options.bsum_weight/2*G.fac{m}{k};
                                    B{m}{k} = B{m}{k} + options.bsum_weight/2*eye(size(B{m}{k}));
                                end
                                last_m(p) = 2;
                                if Z.constrained_modes(m) && iter >= options.iter_start_PAR2Bkconstraint
                                    B{m}{k} = B{m}{k} + rho{m}(k)/2*eye(size(B{m}{k}));
                                end
                                L{m}{k} = chol(B{m}{k},'lower'); %precompute Cholesky decomposition 
                            end
                            [inner_iters,G] = ADMM_B_Parafac2(Z,G,iter,A{m},L{m},m,p,rho{m},options);
                            out.innerIters(m,iter)= inner_iters;
                            for k=1:length(Z.size{Z.modes{p}(2)})
                                G_transp_G{m}{k} = G.fac{m}{k}'*G.fac{m}{k}; % update G transposed G for mth mode
                            end
                        else % third Parafac2 mode, can be coupled
                            for k=1:length(Z.size{Z.modes{p}(2)})  
                                A{m}{k} = Z.weights(p)*diag(G.fac{Z.modes{p}(1)}'*Z.object{p}{k}*G.fac{Z.modes{p}(2)}{k});
                                C{m}{k} = G_transp_G{Z.modes{p}(1)}.*G_transp_G{Z.modes{p}(2)}{k};
                                rho{m}(k) = trace(C{m}{k})/size(C{m}{k},1);
                                B{m}{k} = Z.weights(p)* C{m}{k};
                                if isfield(Z,'ridge')
                                    B{m}{k} = B{m}{k} + Z.ridge(m)*eye(size(B{m}{k}));
                                end
                                %last_mttkrp{p}{k} = A{m}{k}*1/Z.weights(p);
                                last_m(p) = 3;
                                if options.bsum
                                    A{m}{k} = A{m}{k} + options.bsum_weight/2*G.fac{m}(k,:)';
                                    B{m}{k} = B{m}{k} + options.bsum_weight/2*eye(size(B{m}{k}));
                                end
                                if coupl_id==0 %mode is not coupled
                                    if (Z.constrained_modes(m)==0) % mode is not constrained
                                        G.fac{m}(k,:) = (B{m}{k}\A{m}{k})'; % ALS update (row-wise)
                                        inner_iters = 1;
                                    else  % is constrained, use ADMM (here only precomputations)
                                        B{m}{k} = B{m}{k} + rho{m}(k)/2*eye(size(B{m}{k}));
                                        L{m}{k} = chol(B{m}{k}','lower'); %precompute Cholesky decomposition
                                    end
                                end
                            end
                            if Z.constrained_modes(m) && coupl_id == 0 % is constrained and not coupled, use ADMM
                                [inner_iters,lbfgsb_iterations,G] = ADMM_constrained_only(Z,G,nb_modes,A{m},L{m},m,p,rho,options);
                            end
                            out.innerIters(m,iter)= inner_iters;
                        end
                    end
                end
            end

                if coupl_id~=0 %modes are coupled: use "coupled ADMM"
                    ctype = Z.coupling.coupling_type(coupl_id); %type of linear coupling
                    switch ctype
                        case 0 %exact coupling
                            for m=coupled_modes
                                p = which_p(m);
                                if strcmp(Z.loss_function{p},'Frobenius')
                                    if strcmp(Z.model{p},'PAR2') && 3 == find(Z.modes{p}==m) % third Parafac2 mode
                                        for k=1:length(Z.size{Z.modes{p}(2)}) 
                                            B{m}{k} = B{m}{k} + rho{m}(k)/2*eye(size(B{m}{k})); %for the coupling
                                            if Z.constrained_modes(m) %mode is constrained
                                               B{m}{k} = B{m}{k} + rho{m}(k)/2*eye(size(B{m}{k})); %for the constraint
                                            end
                                            L{m}{k} = chol(B{m}{k}','lower'); %precompute Cholesky decomposition
                                        end
                                    else
                                        B{m} = B{m} + rho{m}/2* eye(size(B{m})); % for the coupling
                                        if Z.constrained_modes(m) %mode is constrained
                                            B{m} = B{m} + rho{m}/2*eye(size(B{m}));
                                        end
                                        L{m} = chol(B{m}','lower'); %precompute Cholesky decomposition of B (only works in the chase when rho does not change between inner iterations)
                                    end
                                end
                            end
                            [inner_iters,lbfgsb_iterations,G] = ADMM_coupled_case0(Z,G,nb_modes,which_p,m,lscalar,uscalar,fh,gh,A,L,coupled_modes,coupl_id,rho,options);
                        case 1 % mode is linear coupled with trafo matrix from left, use ADMM with silvester equation
                            for m=coupled_modes
                                p = which_p(m);
                                if strcmp(Z.loss_function{p},'Frobenius')
                                    B2{m} = rho{m}/2* Z.coupling.coupl_trafo_matrices{m}'*Z.coupling.coupl_trafo_matrices{m}; % precompute????
                                    if Z.constrained_modes(m) %mode is constrained 
                                        B2{m} = B2{m} + rho{m}/2*eye(size(B2{m}));
                                    end
                                end
                            end
                            [inner_iters,lbfgsb_iterations,G] = ADMM_coupled_case1(Z,G,nb_modes,which_p,m,A,B,B2,coupled_modes,coupl_id,rho,options);
                        case 2
                            for m=coupled_modes
                                p = which_p(m);
                                if strcmp(Z.loss_function{p},'Frobenius')
                                    B{m} = B{m} + rho{m}/2* Z.coupling.coupl_trafo_matrices{m}*Z.coupling.coupl_trafo_matrices{m}'; % precompute????; % for the coupling
                                    if Z.constrained_modes(m) %mode is constrained
                                        B{m} = B{m} + rho{m}/2*eye(size(B{m}));
                                    end
                                    L{m} = chol(B{m}','lower'); %precompute Cholesky decomposition of B (only works in the chase when rho does not change between inner iterations)
                                end
                            end
                            [inner_iters,lbfgsb_iterations,G] = ADMM_coupled_case2(Z,G,nb_modes,which_p,m,A,L,coupled_modes,coupl_id,rho,options);
                        case 3
                            for m=coupled_modes
                                p = which_p(m);
                                if strcmp(Z.loss_function{p},'Frobenius')
                                    B{m} = B{m} + rho{m}/2* eye(size(B{m})); % for the coupling
                                    if Z.constrained_modes(m) %mode is constrained
                                        B{m} = B{m} + rho{m}/2*eye(size(B{m}));
                                    end
                                    L{m} = chol(B{m}','lower'); %precompute Cholesky decomposition of B (only works in the chase when rho does not change between inner iterations)
                                end
                            end
                            [inner_iters,lbfgsb_iterations,G] = ADMM_coupled_case3(Z,G,nb_modes,which_p,m,A,L,coupled_modes,coupl_id,rho,options);
                        case 4
                            for m=coupled_modes
                                p = which_p(m);
                                if strcmp(Z.loss_function{p},'Frobenius')
                                    B{m} = B{m} + rho{m}/2* eye(size(B{m})); % for the coupling
                                    if Z.constrained_modes(m) %mode is constrained
                                        B{m} = B{m} + rho{m}/2*eye(size(B{m}));
                                    end
                                    L{m} = chol(B{m}','lower'); %precompute Cholesky decomposition of B (only works in the chase when rho does not change between inner iterations)
                                end
                            end
                            [inner_iters,lbfgsb_iterations,G] = ADMM_coupled_case4(Z,G,nb_modes,which_p,m,A,L,coupled_modes,coupl_id,rho,options);
                    end
                    out.innerIters(coupled_modes,iter)= inner_iters;
                    for m=coupled_modes
                        p = which_p(m);
                        if strcmp(Z.loss_function{p},'Frobenius')
                            G_transp_G{m} = G.fac{m}'*G.fac{m}; % update G transposed G for mth mode
                        else
                            out.lbfgsb_iterations{m,iter} = lbfgsb_iterations{m};
                            for r=1:size(G.fac{m},2)
                                sum_column_norms_sqr(m,1) = norm(G.fac{m}(:,r))^2;
                            end
                        end
                    end
                end
            %end
        end
        
        
       f_tensors_old = f_tensors;
       f_couplings_old = f_couplings;
       f_constraints_old = f_constraints;
       f_PAR2_couplings_old = f_PAR2_couplings;
       [f_tensors,f_couplings,f_constraints,f_PAR2_couplings] = CMTF_AOADMM_func_eval(Z,G,G_transp_G,Znorm_const,last_mttkrp,last_had,last_m,fh);
       f_total = f_tensors+f_couplings+f_constraints + f_PAR2_couplings;
       func_val(iter+1) = f_tensors;
       func_coupl(iter+1) = f_couplings;
       func_constr(iter+1) = f_constraints;
       func_PAR2_coupl(iter+1) = f_PAR2_couplings;
       time_at_it(iter+1) = toc(tstart);
       stop = cmtf.utils.evaluate_stopping_conditions(f_tensors,f_couplings,f_constraints,f_PAR2_couplings,f_tensors_old,f_couplings_old,f_constraints_old,f_PAR2_couplings_old,options);

        %display
        if strcmp(options.Display,'iter') && mod(iter,options.DisplayIters)==0
            fprintf(1,'%6d %12f %12f %12f %17f %12f\n', iter, f_total, f_tensors, f_couplings,f_constraints,f_PAR2_couplings);
        end
        iter = iter+1;
        
        % catch
        %     illconditioned = 1;
        %     fprintf('Stopped due to illconditioned linear system')
        %     break
        % end
    end
    % which condition caused stop?
    exit_flag = cmtf.utils.make_exit_flag(iter,f_tensors,f_couplings,f_constraints,f_PAR2_couplings,options,illconditioned);
    %save output
    out.f_tensors = f_tensors;
    out.f_couplings = f_couplings;
    out.f_constraints = f_constraints;
    out.f_PAR2_couplings = f_PAR2_couplings;
    out.exit_flag = exit_flag;
    out.OuterIterations = iter-1;
    out.func_val_conv = func_val;
    out.func_coupl_conv = func_coupl;
    out.func_constr_conv = func_constr;
    out.func_PAR2_coupl = func_PAR2_coupl;
    out.time_at_it = time_at_it;


    %display final
    if strcmp(options.Display,'iter') || strcmp(options.Display,'final')
        fprintf(1,'%6d %12f %12f %12f %12f %12f\n', iter-1, f_total, f_tensors, f_couplings,f_constraints,f_PAR2_couplings);
    end
    
    warning(s); 

end

