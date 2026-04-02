% Test missing data handling via EM imputation for a single CP tensor.
%
% Creates a noiseless rank-3 CP tensor of size [20 30 40], masks ~20% of
% entries at random, and verifies that the AO-ADMM EM solver recovers the
% factors and imputed entries to high accuracy.

function tests = test_missing_data_CP
    tests = functiontests(localfunctions);
end

function test_CP_missing_data_recovery(testCase)

    rng(42, 'twister');

    % Problem dimensions
    sz     = {20, 30, 40};
    modes  = {[1 2 3]};
    R      = 3;

    % Single uncoupled CP tensor
    model{1}              = 'CP';
    loss_function{1}      = 'Frobenius';
    loss_function_param{1} = [];
    coupling.lin_coupled_modes  = [0 0 0];
    coupling.coupling_type      = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    constrained_modes = [0 0 0];
    constraints       = cell(3,1);
    weights           = [1];

    % Build Z
    Z.loss_function       = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model               = model;
    Z.modes               = modes;
    Z.size                = sz;
    Z.coupling            = coupling;
    Z.constrained_modes   = constrained_modes;
    Z.constraints         = constraints;
    Z.weights             = weights;

    % Create noiseless synthetic data
    lambdas_data   = {ones(1, R)};
    noise          = 0;
    normalize_cols = 0;
    distr_data     = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    [X, Atrue] = cmtf.utils.create_coupled_data('model', model, 'size', sz, ...
        'modes', modes, 'lambdas', lambdas_data, 'noise', noise, ...
        'coupling', coupling, 'normalize_columns', normalize_cols, ...
        'distr_data', distr_data, 'loss_function', Z.loss_function);

    % Normalize data
    normZ      = norm(X{1});
    Z.object{1} = X{1} / normZ;

    % Build missing mask: ~20% missing
    full_size = sz{1} * sz{2} * sz{3};
    miss_mask       = true(sz{1}, sz{2}, sz{3});
    miss_mask(randperm(full_size, round(0.2 * full_size))) = false;
    Z.miss{1} = miss_mask;

    % Initialization
    init_options.lambdas_init = {ones(1, R)};
    init_options.nvecs        = 0;
    init_options.distr        = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    init_options.normalize    = 1;
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

    % Solver options
    options.Display              = 'no';
    options.DisplayIters         = 10;
    options.MaxOuterIters        = 4000;
    options.MaxInnerIters        = 5;
    options.AbsFuncTol           = 1e-7;
    options.OuterRelTol          = 1e-8;
    options.innerRelPrTol_coupl  = 1e-5;
    options.innerRelDualTol_coupl = 1e-5;
    options.innerRelPrTol_constr = 1e-5;
    options.innerRelDualTol_constr = 1e-5;
    options.bsum                 = 0;
    options.eps_log              = 1e-10;

    % Run solver
    [Zhat, ~, ~, out] = cmtf.aoadmm.cmtf_AOADMM(Z, 'alg_options', options, ...
        'init', init_fac, 'init_options', init_options);

    % Tensor fit on observed entries (Z.object{1} already has missing = 0)
    obs_mask = Z.miss{1};
    X_obs    = double(Z.object{1});
    X_rec    = double(full(Zhat{1}));
    fit_obs  = 100 * (1 - norm(X_obs(obs_mask) - X_rec(obs_mask))^2 / norm(X_obs(obs_mask))^2);

    % Factor match score
    true_ktensor = ktensor(lambdas_data{1}' ./ normZ, Atrue(modes{1}));
    FMS = score(Zhat{1}, true_ktensor, 'lambda_penalty', false);

    % Assertions
    testCase.verifyEqual(out.f_tensors, 0, 'AbsTol', 1e-3);
    testCase.verifyEqual(fit_obs, 100, 'AbsTol', 1e-1);
    testCase.verifyEqual(FMS, 1, 'AbsTol', 1e-2);

    % EM-specific output
    testCase.verifyTrue(isfield(out, 'f_rel_missing_final'), ...
        'out.f_rel_missing_final should be present when Z.miss is provided');
    testCase.verifyLessThan(out.f_rel_missing_final, options.OuterRelTol);
end

function test_fully_observed_mask_matches_no_mask(testCase)
    % When Z.miss is all-true (no entries missing), result must be the
    % same as running without Z.miss.

    rng(7, 'twister');

    sz     = {15, 20, 25};
    modes  = {[1 2 3]};
    R      = 2;

    model{1}              = 'CP';
    loss_function{1}      = 'Frobenius';
    loss_function_param{1} = [];
    coupling.lin_coupled_modes   = [0 0 0];
    coupling.coupling_type       = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    constrained_modes = [0 0 0];
    constraints       = cell(3,1);
    weights           = [1];

    Z_base.loss_function       = loss_function;
    Z_base.loss_function_param = loss_function_param;
    Z_base.model               = model;
    Z_base.modes               = modes;
    Z_base.size                = sz;
    Z_base.coupling            = coupling;
    Z_base.constrained_modes   = constrained_modes;
    Z_base.constraints         = constraints;
    Z_base.weights             = weights;

    lambdas_data = {ones(1, R)};
    distr_data   = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    [X] = cmtf.utils.create_coupled_data('model', model, 'size', sz, ...
        'modes', modes, 'lambdas', lambdas_data, 'noise', 0, ...
        'coupling', coupling, 'normalize_columns', 0, ...
        'distr_data', distr_data, 'loss_function', Z_base.loss_function);

    normZ = norm(X{1});
    Z_base.object{1} = X{1} / normZ;

    init_options.lambdas_init = {ones(1, R)};
    init_options.nvecs        = 0;
    init_options.distr        = distr_data;
    init_options.normalize    = 1;
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z_base, 'init_options', init_options);

    options.Display              = 'no';
    options.DisplayIters         = 10;
    options.MaxOuterIters        = 500;
    options.MaxInnerIters        = 5;
    options.AbsFuncTol           = 1e-7;
    options.OuterRelTol          = 1e-6;
    options.innerRelPrTol_coupl  = 1e-5;
    options.innerRelDualTol_coupl = 1e-5;
    options.innerRelPrTol_constr = 1e-5;
    options.innerRelDualTol_constr = 1e-5;
    options.bsum                 = 0;
    options.eps_log              = 1e-10;

    % Baseline: no Z.miss
    [~, ~, ~, out_base] = cmtf.aoadmm.cmtf_AOADMM(Z_base, 'alg_options', options, ...
        'init', init_fac, 'init_options', init_options);

    % With all-observed mask: should not change behaviour
    Z_mask = Z_base;
    Z_mask.object{1} = Z_base.object{1};   % reset (cmtf_AOADMM modifies Z internally)
    Z_mask.miss{1} = true(sz{1}, sz{2}, sz{3});

    [~, ~, ~, out_mask] = cmtf.aoadmm.cmtf_AOADMM(Z_mask, 'alg_options', options, ...
        'init', init_fac, 'init_options', init_options);

    testCase.verifyEqual(out_mask.f_tensors, out_base.f_tensors, 'AbsTol', 1e-6, ...
        'All-observed mask should give same f_tensors as no-mask run');
    testCase.verifyFalse(isfield(out_base, 'f_rel_missing_final'), ...
        'f_rel_missing_final should NOT be in out when Z.miss is absent');
    testCase.verifyTrue(isfield(out_mask, 'f_rel_missing_final'), ...
        'f_rel_missing_final should be in out when Z.miss is present');
end
