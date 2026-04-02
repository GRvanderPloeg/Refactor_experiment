% Test missing data handling via EM imputation for a single PAR2 tensor.
%
% Creates a noiseless rank-3 PARAFAC2 tensor (20 slices of size [20 x 30]),
% masks ~20% of entries per slice, and verifies that the AO-ADMM EM solver
% recovers the shared factor A and the weight matrix C to high accuracy.

function tests = test_missing_data_PAR2
    tests = functiontests(localfunctions);
end

function test_PAR2_missing_data_recovery(testCase)

    rng(55, 'twister');

    K  = 20;   % number of slices
    I  = 20;   % rows (shared mode)
    Jk = 30;   % columns (varying mode, uniform here)
    R  = 3;

    sz    = {I, Jk * ones(1, K), K};
    modes = {[1 2 3]};

    model{1}              = 'PAR2';
    loss_function{1}      = 'Frobenius';
    loss_function_param{1} = [];
    coupling.lin_coupled_modes   = [0 0 0];
    coupling.coupling_type       = [];
    coupling.coupl_trafo_matrices = cell(3,1);
    constrained_modes = [0 0 0];
    constraints       = cell(3,1);
    weights           = [1];

    Z.loss_function       = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model               = model;
    Z.modes               = modes;
    Z.size                = sz;
    Z.coupling            = coupling;
    Z.constrained_modes   = constrained_modes;
    Z.constraints         = constraints;
    Z.weights             = weights;

    lambdas_data = {ones(1, R)};
    distr_data   = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) rand(x,y)+0.1};
    [X, Atrue] = cmtf.utils.create_coupled_data('model', model, 'size', sz, ...
        'modes', modes, 'lambdas', lambdas_data, 'noise', 0, ...
        'coupling', coupling, 'normalize_columns', 0, ...
        'distr_data', distr_data, 'loss_function', Z.loss_function);

    % Normalize
    normZ = 0;
    for k = 1:K
        normZ = normZ + norm(X{1}{k}, 'fro')^2;
    end
    normZ = sqrt(normZ);
    Z.object{1} = X{1};
    for k = 1:K
        Z.object{1}{k} = X{1}{k} / normZ;
    end

    % Build per-slice missing masks: ~20% missing per slice
    Z.miss{1} = cell(K, 1);
    for k = 1:K
        n_entries  = I * Jk;
        mask_k     = true(I, Jk);
        mask_k(randperm(n_entries, round(0.2 * n_entries))) = false;
        Z.miss{1}{k} = mask_k;
    end

    init_options.lambdas_init = {ones(1, R)};
    init_options.nvecs        = 0;
    init_options.distr        = distr_data;
    init_options.normalize    = 1;
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

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

    [Zhat, ~, ~, out] = cmtf.aoadmm.cmtf_AOADMM(Z, 'alg_options', options, ...
        'init', init_fac, 'init_options', init_options);

    % Factor match scores for A and C
    FMS_A = score(ktensor(ones(R,1), Zhat{1}.A), ...
                  ktensor(ones(R,1), Atrue{1}), 'lambda_penalty', false);
    FMS_C = score(ktensor(ones(R,1), Zhat{1}.C), ...
                  ktensor(ones(R,1), Atrue{3}), 'lambda_penalty', false);

    % Assertions
    testCase.verifyEqual(out.f_tensors, 0, 'AbsTol', 1e-3);
    testCase.verifyEqual(FMS_A, 1, 'AbsTol', 1e-2);
    testCase.verifyEqual(FMS_C, 1, 'AbsTol', 1e-2);

    % EM-specific output
    testCase.verifyTrue(isfield(out, 'f_rel_missing_final'), ...
        'out.f_rel_missing_final should be present for PAR2 missing data');
    testCase.verifyLessThan(out.f_rel_missing_final, options.OuterRelTol);
end
