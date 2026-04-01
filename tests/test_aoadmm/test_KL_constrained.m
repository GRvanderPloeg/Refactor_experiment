% Test that compute_gen_f_g handles constrained=true (lines 11-12).
%
% Two coupled CP tensors with KL loss and coupling type 0.  The two coupled
% modes are also constrained to be non-negative.  This forces execution
% through ADMM_coupled_case0 -> lbfgsb_update(..., constrained=1, coupling_type=0, ...)
% -> compute_gen_f_g with constrained==true, covering lines 11-12.

function tests = test_KL_constrained
    tests = functiontests(localfunctions);
end

function test_workflow(testCase)

    rng(42, 'twister');

    sz     = {20, 15, 12, 20, 10};
    P      = 2;
    lambdas_data = {[1 1], [1 1]};
    modes  = {[1 2 3], [4 5]};
    noise  = 0.1;
    shape  = 1; scale = 1;
    distr_data = {@(x,y) gamrnd(shape,scale,x,y), @(x,y) gamrnd(shape,scale,x,y), ...
                  @(x,y) gamrnd(shape,scale,x,y), @(x,y) gamrnd(shape,scale,x,y), ...
                  @(x,y) gamrnd(shape,scale,x,y)};
    normalize_columns = 0;

    model{1} = 'CP';
    model{2} = 'CP';

    % Exact coupling (type 0) on mode 1 and 4
    coupling.lin_coupled_modes    = [1 0 0 1 0];
    coupling.coupling_type        = [0];
    coupling.coupl_trafo_matrices = cell(5, 1);

    loss_function{1}       = 'KL';
    loss_function{2}       = 'KL';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    cmtf.utils.check_data_input(sz, modes, lambdas_data, coupling, loss_function, model);

    init_options.lambdas_init = {[1 1], [1 1]};
    init_options.nvecs        = 0;
    init_options.distr        = {@(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y), ...
                                  @(x,y) rand(x,y), @(x,y) rand(x,y)};
    init_options.normalize    = 1;

    % Non-negativity constraint on the two coupled modes
    constrained_modes = [1 0 0 1 0];
    constraints       = cell(length(constrained_modes), 1);
    constraints{1}    = {'non-negativity'};
    constraints{4}    = {'non-negativity'};

    weights = [1/2 1/2];

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
    options.MaxOuterIters          = 200;
    options.MaxInnerIters          = 3;
    options.AbsFuncTol             = 1e-4;
    options.OuterRelTol            = 1e-6;
    options.innerRelPrTol_coupl    = 1e-3;
    options.innerRelPrTol_constr   = 1e-3;
    options.innerRelDualTol_coupl  = 1e-3;
    options.innerRelDualTol_constr = 1e-3;
    options.bsum                   = 0;
    options.eps_log                = 1e-10;
    options.lbfgsb_options         = lbfgsb_options;

    [~, ~, ~, out] = cmtf.aoadmm.cmtf_AOADMM(Z, ...
        'alg_options', options, 'init', init_fac, 'init_options', init_options);

    verifyTrue(testCase, isfinite(out.f_tensors));
end
