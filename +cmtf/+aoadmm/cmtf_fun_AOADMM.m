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
        fprintf(1,' Iter  f total      f tensors      f couplings    f constraints    f PAR2 couplings\n');
        fprintf(1,'------ ------------ -------------  -------------- ---------------- ----------------\n');
    end

    if strcmp(options.Display,'iter')
        fprintf(1,'%6d %12f %12f %12f %17f %12f\n', 0, f_total, f_tensors, f_couplings,f_constraints, f_PAR2_couplings);
    end
    iter = 1;
    
    [G_transp_G,sum_column_norms_sqr] = init_G_transp_G(Z,G,which_p,nb_modes);

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
                        [A{m},B{m},C{m},rho{m},last_mttkrp,last_had,last_m] = ...
                            precompute_mode_cp(Z,G,G_transp_G,sum_column_norms_sqr,...
                                               m,p,last_mttkrp,last_had,last_m,options);
                        if coupl_id==0 %modes are not coupled
                            [G,G_transp_G,sum_column_norms_sqr,inner_iters,lbfgsb_iterations] = ...
                                dispatch_uncoupled(Z,G,nb_modes,G_transp_G,sum_column_norms_sqr,...
                                                   A{m},B{m},[],rho{m},m,p,rho,...
                                                   lscalar,uscalar,fh,gh,options);
                            out.innerIters(m,iter) = inner_iters;
                            if ~strcmp(Z.loss_function{p},'Frobenius')
                                out.lbfgsb_iterations{m,iter} = lbfgsb_iterations;
                            end
                        end
                    end
                elseif strcmp(Z.model{p},'PAR2')
                    for m=coupled_modes(which_p(coupled_modes)==p)
                        [A{m},B{m},C{m},rho{m},L{m},G,par2_iters,...
                         last_mttkrp,last_had,last_m] = ...
                            precompute_mode_par2(Z,G,G_transp_G,m,p,iter,coupl_id,...
                                                 last_mttkrp,last_had,last_m,options);
                        mode_pos = find(Z.modes{p}==m);
                        if mode_pos == 2 % second PAR2 mode (the "funny" mode)
                            [inner_iters,G] = ADMM_B_Parafac2(Z,G,iter,A{m},L{m},m,p,rho{m},options);
                            out.innerIters(m,iter) = inner_iters;
                            for k=1:length(Z.size{Z.modes{p}(2)})
                                G_transp_G{m}{k} = G.fac{m}{k}'*G.fac{m}{k};
                            end
                        else % mode_pos 1 or 3
                            if coupl_id==0
                                if ~isempty(par2_iters) % mode 3 unconstrained: update done in precompute
                                    inner_iters = par2_iters;
                                else
                                    [G,G_transp_G,sum_column_norms_sqr,inner_iters,~] = ...
                                        dispatch_uncoupled(Z,G,nb_modes,G_transp_G,...
                                                           sum_column_norms_sqr,...
                                                           A{m},B{m},L{m},rho{m},...
                                                           m,p,rho,lscalar,uscalar,fh,gh,options);
                                end
                            end
                            if mode_pos == 1
                                G_transp_G{m} = G.fac{m}'*G.fac{m};
                            end
                            out.innerIters(m,iter) = inner_iters;
                        end
                    end
                end
            end

                if coupl_id~=0 %modes are coupled: use "coupled ADMM"
                    ctype = Z.coupling.coupling_type(coupl_id); %type of linear coupling
                    [inner_iters,lbfgsb_iterations,G] = ADMM_coupled(Z,G,nb_modes,which_p,ctype,lscalar,uscalar,fh,gh,A,B,coupled_modes,coupl_id,rho,options);
                    out.innerIters(coupled_modes,iter) = inner_iters;
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

