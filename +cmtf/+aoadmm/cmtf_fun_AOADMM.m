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
    has_missing = isfield(Z, 'miss') && any(~cellfun(@isempty, Z.miss));
    f_rel_missing = Inf;
    G_transp_G = cell(nb_modes,1);
    A = cell(nb_modes,1);
    C = cell(nb_modes,1);
    B = cell(nb_modes,1);
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
        if has_missing
            fprintf(1,' Iter  f total      f tensors      f couplings    f constraints    f PAR2 couplings  f_rel_miss\n');
        else
            fprintf(1,' Iter  f total      f tensors      f couplings    f constraints    f PAR2 couplings\n');
        end
        fprintf(1,'------ ------------ -------------  -------------- ---------------- ----------------\n');
    end

    if strcmp(options.Display,'iter')
        if has_missing
            fprintf(1,'%6d %12f %12f %12f %17f %12f %12f\n', 0, f_total, f_tensors, f_couplings, f_constraints, f_PAR2_couplings, f_rel_missing);
        else
            fprintf(1,'%6d %12f %12f %12f %17f %12f\n', 0, f_total, f_tensors, f_couplings,f_constraints, f_PAR2_couplings);
        end
    end
    iter = 1;
    
    [G_transp_G,sum_column_norms_sqr] = init_G_transp_G(Z,G,which_p,nb_modes);

    s = warning('error', 'MATLAB:nearlySingularMatrix');
    illconditioned = 0;
    
    stop = false;
    while(iter<=options.MaxOuterIters && ~stop)
        %try

        % --- Alternating Optimization sweep ---
        % Iterate over all coupling groups. coupl_id=0 collects all uncoupled
        % modes so they are processed with the same loop structure.
        for coupl_id=couplings
            coupled_modes = find(Z.coupling.lin_coupled_modes==coupl_id); % modes belonging to this coupling group

            % Loop over the tensors that contain at least one mode in this
            % coupling group. Each tensor is processed independently here;
            % the inter-tensor coupling constraint is enforced afterwards.
            for p=unique(which_p(coupled_modes))
                if strcmp(Z.model{p},'CP')

                    % For each mode of tensor p that belongs to this coupling group,
                    % precompute the ADMM subproblem ingredients (normal equations
                    % LHS A, RHS B, Gram matrix C, and penalty parameter rho).
                    % Modes within the same tensor must be updated sequentially
                    % because each update changes G, which invalidates cached MTTKRP.
                    for m=coupled_modes(which_p(coupled_modes)==p)
                        [A{m},B{m},C{m},rho{m},last_mttkrp,last_had,last_m] = ...
                            precompute_mode_cp(Z,G,G_transp_G,sum_column_norms_sqr,...
                                               m,p,last_mttkrp,last_had,last_m,options);

                        if coupl_id==0 % uncoupled mode: update factor immediately
                            % dispatch_uncoupled solves the ADMM subproblem for
                            % mode m and updates G, G_transp_G, and sum_column_norms_sqr.
                            [G,G_transp_G,sum_column_norms_sqr,inner_iters,lbfgsb_iterations] = ...
                                dispatch_uncoupled(Z,G,nb_modes,G_transp_G,sum_column_norms_sqr,...
                                                   A{m},B{m},[],rho{m},m,p,rho,...
                                                   lscalar,uscalar,fh,gh,options);
                            out.innerIters(m,iter) = inner_iters;
                            if ~strcmp(Z.loss_function{p},'Frobenius')
                                out.lbfgsb_iterations{m,iter} = lbfgsb_iterations;
                            end
                        end
                        % If coupl_id~=0, the factor update is deferred to the
                        % coupled ADMM step below, once all tensors are precomputed.
                    end

                elseif strcmp(Z.model{p},'PAR2')

                    % PARAFAC2 modes are handled per position within the tensor.
                    % precompute_mode_par2 returns the orthogonal projection
                    % matrices L{m} (one per slice) and may directly update G
                    % for the unconstrained mode (mode position 3).
                    for m=coupled_modes(which_p(coupled_modes)==p)
                        [A{m},B{m},C{m},rho{m},L{m},G,par2_iters,...
                         last_mttkrp,last_had,last_m] = ...
                            precompute_mode_par2(Z,G,G_transp_G,m,p,iter,coupl_id,...
                                                 last_mttkrp,last_had,last_m,options);
                        mode_pos = find(Z.modes{p}==m); % position of m within tensor p (1, 2, or 3)

                        if mode_pos == 2
                            % Mode 2 is the varying-size PARAFAC2 mode (the "funny" mode).
                            % Its factor matrices B_k differ per slice, so a dedicated
                            % ADMM solver is required. G_transp_G must be recomputed
                            % slice-by-slice after the update.
                            [inner_iters,G] = ADMM_B_Parafac2(Z,G,iter,A{m},L{m},m,p,rho{m},options);
                            out.innerIters(m,iter) = inner_iters;
                            for k=1:length(Z.size{Z.modes{p}(2)})
                                G_transp_G{m}{k} = G.fac{m}{k}'*G.fac{m}{k};
                            end

                        else % mode_pos 1 (shared factor) or mode_pos 3 (e.g. time mode)
                            if coupl_id==0 % uncoupled: update factor now
                                if ~isempty(par2_iters)
                                    % Mode 3 unconstrained: precompute_mode_par2 already
                                    % performed a closed-form update; just record iterations.
                                    inner_iters = par2_iters;
                                else
                                    % Mode 1 or constrained mode 3: solve via ADMM.
                                    [G,G_transp_G,sum_column_norms_sqr,inner_iters,~] = ...
                                        dispatch_uncoupled(Z,G,nb_modes,G_transp_G,...
                                                           sum_column_norms_sqr,...
                                                           A{m},B{m},L{m},rho{m},...
                                                           m,p,rho,lscalar,uscalar,fh,gh,options);
                                end
                            end
                            if mode_pos == 1
                                % Mode 1 is a full factor matrix; refresh its Gram matrix.
                                G_transp_G{m} = G.fac{m}'*G.fac{m};
                            end
                            out.innerIters(m,iter) = inner_iters;
                        end
                    end
                end
            end % loop over tensors p

            if coupl_id~=0
                % Coupled modes: enforce the linear coupling constraint across tensors
                % by running a joint ADMM step. A and B precomputed above enter as
                % the per-mode subproblem data; the coupling type determines the
                % constraint structure (e.g. equality, linear transform).
                ctype = Z.coupling.coupling_type(coupl_id);
                [inner_iters,lbfgsb_iterations,G] = ADMM_coupled(Z,G,nb_modes,which_p,ctype,lscalar,uscalar,fh,gh,A,B,coupled_modes,coupl_id,rho,options);
                out.innerIters(coupled_modes,iter) = inner_iters;

                % After the coupled update, refresh the cached quantities used in
                % subsequent MTTKRP and objective evaluations.
                for m=coupled_modes
                    p = which_p(m);
                    if strcmp(Z.loss_function{p},'Frobenius')
                        G_transp_G{m} = G.fac{m}'*G.fac{m};
                    else
                        % For non-Frobenius losses G_transp_G is not needed, but
                        % sum_column_norms_sqr is used in the L-BFGS-B line search.
                        out.lbfgsb_iterations{m,iter} = lbfgsb_iterations{m};
                        for r=1:size(G.fac{m},2)
                            sum_column_norms_sqr(m,1) = norm(G.fac{m}(:,r))^2;
                        end
                    end
                end
            end
            %end
        end % loop over coupling groups

        % EM imputation: update missing entries with the current low-rank model.
        % After this step, Z.object{p}(missing) == model(missing), so the
        % Frobenius residual at missing positions is exactly zero and the
        % objective measures fit on observed entries only.
        if has_missing
            num_sq = 0;
            den_sq = 0;
            for p = 1:P
                if isempty(Z.miss{p}), continue; end
                if strcmp(Z.model{p}, 'CP')
                    M_full = double(full(ktensor(G.fac(Z.modes{p}))));
                    miss_mask = ~Z.miss{p};
                    % tensor objects do not support 3-D logical indexing;
                    % extract as double, modify, then wrap back.
                    if isa(Z.object{p}, 'tensor')
                        tmp = double(Z.object{p});
                        old_vals = tmp(miss_mask);
                        new_vals = M_full(miss_mask);
                        tmp(miss_mask) = new_vals;
                        Z.object{p} = tensor(tmp);
                        Znorm_const{p} = norm(Z.object{p})^2;
                    else
                        old_vals = Z.object{p}(miss_mask);
                        new_vals = M_full(miss_mask);
                        Z.object{p}(miss_mask) = new_vals;
                        Znorm_const{p} = norm(Z.object{p}(:))^2;
                    end
                    num_sq = num_sq + sum((new_vals - old_vals).^2);
                    den_sq = den_sq + sum(old_vals.^2);
                elseif strcmp(Z.model{p}, 'PAR2')
                    m1 = Z.modes{p}(1); m2 = Z.modes{p}(2); m3 = Z.modes{p}(3);
                    Znorm_const{p} = 0;
                    for k = 1:length(Z.object{p})
                        M_k = G.fac{m1} * diag(G.fac{m3}(k,:)) * G.fac{m2}{k}';
                        miss_k = ~Z.miss{p}{k};
                        old_k = Z.object{p}{k}(miss_k);
                        new_k = M_k(miss_k);
                        num_sq = num_sq + sum((new_k - old_k).^2);
                        den_sq = den_sq + sum(old_k.^2);
                        Z.object{p}{k}(miss_k) = new_k;
                        Znorm_const{p} = Znorm_const{p} + norm(Z.object{p}{k}, 'fro')^2;
                    end
                end
                % Invalidate the per-tensor MTTKRP cache so func_eval uses the
                % slow path (recomputes from the updated Z.object{p}).
                last_mttkrp{p} = [];
                last_had{p} = [];
            end
            if den_sq > 0
                f_rel_missing = sqrt(num_sq / den_sq);
            else
                f_rel_missing = sqrt(num_sq);
            end
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
       if has_missing
           stop = stop && (f_rel_missing < options.OuterRelTol);
       end

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
    if has_missing
        out.f_rel_missing_final = f_rel_missing;
    end
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

