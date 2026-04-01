% Test that compute_gen_f_g handles the bsum branch (lines 37-38).
%
% Single uncoupled, unconstrained CP tensor with KL loss and options.bsum=1.
% This routes through the direct lbfgsb path in cmtf_fun_AOADMM (line 151),
% which calls compute_gen_f_g.  With options.bsum true, lines 37-38 execute.

function tests = test_KL_bsum
    tests = functiontests(localfunctions);
end

function test_workflow(testCase)

    rng(11, 'twister');

    sz     = {15, 10, 12};
    P      = 1;
    lambdas_data = {[1 1]};
    modes  = {[1 2 3]};
    noise  = 0.1;
    shape  = 1; scale = 1;
    distr_data = {@(x,y) gamrnd(shape,scale,x,y), @(x,y) gamrnd(shape,scale,x,y), ...
                  @(x,y) gamrnd(shape,scale,x,y)};
    normalize_columns = 0;

    model{1} = 'CP';

    coupling.lin_coupled_modes    = [0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(3, 1);

    loss_function{1}       = 'KL';
    loss_function_param{1} = [];

    cmtf.utils.check_data_input(sz, modes, lambdas_data, coupling, loss_function, model);

    init_options.lambdas_init = {[1 1]};
    init_options.nvecs        = 0;
    init_options.distr        = {@(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y)};
    init_options.normalize    = 1;

    constrained_modes = [0 0 0];
    constraints       = cell(length(constrained_modes), 1);

    weights = [1];

    Z.loss_function       = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model               = model;
    Z.modes               = modes;
    Z.size                = sz;
    Z.coupling            = coupling;
    Z.constrained_modes   = constrained_modes;
    Z.constraints         = constraints;
    Z.weights             = weights;

    [X, ~, ~, ~] = cmtf.utils.create_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, ...
        'lambdas', lambdas_data, 'noise', noise, ...
        'coupling', coupling, 'normalize_columns', normalize_columns, ...
        'distr_data', distr_data, 'loss_function', Z.loss_function);
    for p = 1:P
        Z.object{p} = X{p};
    end

    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

    lbfgsb_options.m            = 5;
    lbfgsb_options.printEvery   = -1;
    lbfgsb_options.maxIts       = 50;
    lbfgsb_options.maxTotalIts  = 500;
    lbfgsb_options.factr        = 1e-6/eps;
    lbfgsb_options.pgtol        = 1e-4;

    options.Display                = 'no';
    options.DisplayIters           = 10;
    options.MaxOuterIters          = 100;
    options.MaxInnerIters          = 3;
    options.AbsFuncTol             = 1e-4;
    options.OuterRelTol            = 1e-6;
    options.innerRelPrTol_coupl    = 1e-3;
    options.innerRelPrTol_constr   = 1e-3;
    options.innerRelDualTol_coupl  = 1e-3;
    options.innerRelDualTol_constr = 1e-3;
    options.bsum                   = 1;        % enables bsum branch in compute_gen_f_g
    options.bsum_weight            = 0.1;
    options.eps_log                = 1e-10;
    options.lbfgsb_options         = lbfgsb_options;

    [~, ~, ~, out] = cmtf.aoadmm.cmtf_AOADMM(Z, ...
        'alg_options', options, 'init', init_fac, 'init_options', init_options);

    verifyTrue(testCase, isfinite(out.f_tensors));
end
