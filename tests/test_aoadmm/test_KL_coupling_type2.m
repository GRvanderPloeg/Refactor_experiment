% Test coupling type 2 (C*H = Delta) with KL divergence loss.
%
% Same coupling geometry as test_script_coupling_type2 but with KL loss.
% Expected to fail until ADMM_coupled_case2 lbfgsb_update call is fixed.

function tests = test_KL_coupling_type2
    tests = functiontests(localfunctions);
end

function test_workflow(testCase)

    rng(123, 'twister');

    sz           = {50, 30, 40, 50, 70};
    P            = 2;
    lambdas_data = {[1 1 1 1], [1 1 1]};
    modes        = {[1 2 3], [4 5]};
    noise        = 0.2;
    distr_data   = {@(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y), ...
                    @(x,y) rand(x,y), @(x,y) rand(x,y)};
    normalize_columns = 0;

    model{1} = 'CP';
    model{2} = 'CP';

    coupling.lin_coupled_modes    = [1 0 0 1 0];
    coupling.coupling_type        = [2];
    coupling.coupl_trafo_matrices = cell(5, 1);
    coupling.coupl_trafo_matrices{1} = [eye(3, 3); zeros(1, 3)];  % 4x3
    coupling.coupl_trafo_matrices{4} = eye(3, 3);                  % 3x3

    loss_function{1}       = 'KL';
    loss_function{2}       = 'KL';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    cmtf.utils.check_data_input(sz, modes, lambdas_data, coupling, loss_function, model);

    init_options.lambdas_init = lambdas_data;
    init_options.nvecs        = 0;
    init_options.distr        = distr_data;
    init_options.normalize    = 1;

    constrained_modes = [0 0 0 0 0];
    constraints       = cell(length(sz), 1);
    weights           = [0.5 0.5];

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
