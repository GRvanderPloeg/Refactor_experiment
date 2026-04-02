% Test input validation for the missing data feature.
%
% Verifies that cmtf_AOADMM raises informative errors for invalid Z.miss
% inputs and that unsupported configurations (non-Frobenius loss, sptensor)
% are caught before the solver runs.

function tests = test_missing_data_validation
    tests = functiontests(localfunctions);
end

% ----- shared setup helpers ------------------------------------------------

function Z = build_base_Z(sz, model_str, loss_str)
    % Minimal Z struct for a single uncoupled tensor.
    Z.model               = {model_str};
    Z.loss_function       = {loss_str};
    Z.loss_function_param = {[]};
    Z.modes               = {1:length(sz)};
    Z.size                = num2cell(sz);
    Z.coupling.lin_coupled_modes   = zeros(1, length(sz));
    Z.coupling.coupling_type       = [];
    Z.coupling.coupl_trafo_matrices = cell(length(sz), 1);
    Z.constrained_modes   = zeros(1, length(sz));
    Z.constraints         = cell(length(sz), 1);
    Z.weights             = [1];
    % Provide a simple random dense data tensor
    data = randn(sz);
    Z.object{1} = tensor(data / norm(tensor(data)));
end

function options = fast_options()
    options.Display              = 'no';
    options.DisplayIters         = 1;
    options.MaxOuterIters        = 2;
    options.MaxInnerIters        = 2;
    options.AbsFuncTol           = 1e-7;
    options.OuterRelTol          = 1e-8;
    options.innerRelPrTol_coupl  = 1e-5;
    options.innerRelDualTol_coupl = 1e-5;
    options.innerRelPrTol_constr = 1e-5;
    options.innerRelDualTol_constr = 1e-5;
    options.bsum                 = 0;
    options.eps_log              = 1e-10;
end

% ----- error-path tests ----------------------------------------------------

function test_error_non_frobenius_loss(testCase)
    % Non-Frobenius loss + Z.miss must error.
    rng(1, 'twister');
    sz = [5, 6, 7];
    Z  = build_base_Z(sz, 'CP', 'KL');
    % KL requires non-negative integer data
    Z.object{1} = tensor(abs(round(10 * randn(sz))));
    Z.miss{1}   = true(sz);
    Z.miss{1}(1) = false;

    init_options.lambdas_init = {[1 1]};
    init_options.nvecs        = 0;
    init_options.distr        = {@(x,y) rand(x,y)+0.1, @(x,y) rand(x,y)+0.1, @(x,y) rand(x,y)+0.1};
    init_options.normalize    = 1;
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

    testCase.verifyError(...
        @() cmtf.aoadmm.cmtf_AOADMM(Z, 'alg_options', fast_options(), ...
            'init', init_fac, 'init_options', init_options), ...
        'cmtf:missingData:nonFrobenius', ...
        'Non-Frobenius loss with Z.miss should produce an error');
end

function test_error_mask_size_mismatch(testCase)
    % Z.miss{p} of wrong size must error.
    rng(2, 'twister');
    sz = [10, 12, 8];
    Z  = build_base_Z(sz, 'CP', 'Frobenius');
    Z.miss{1} = true(sz + 1);  % deliberate size mismatch

    init_options.lambdas_init = {[1 1]};
    init_options.nvecs        = 0;
    init_options.distr        = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    init_options.normalize    = 1;
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

    testCase.verifyError(...
        @() cmtf.aoadmm.cmtf_AOADMM(Z, 'alg_options', fast_options(), ...
            'init', init_fac, 'init_options', init_options), ...
        'cmtf:missingData:maskSizeMismatch', ...
        'Size-mismatched Z.miss should produce an error');
end

function test_error_non_logical_mask(testCase)
    % Z.miss{p} that is double (not logical) must error.
    rng(3, 'twister');
    sz = [10, 12, 8];
    Z  = build_base_Z(sz, 'CP', 'Frobenius');
    Z.miss{1} = double(true(sz));   % double, not logical

    init_options.lambdas_init = {[1 1]};
    init_options.nvecs        = 0;
    init_options.distr        = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    init_options.normalize    = 1;
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

    testCase.verifyError(...
        @() cmtf.aoadmm.cmtf_AOADMM(Z, 'alg_options', fast_options(), ...
            'init', init_fac, 'init_options', init_options), ...
        'cmtf:missingData:maskNotLogical', ...
        'Non-logical Z.miss should produce an error');
end

function test_error_PAR2_mask_not_cell(testCase)
    % For PAR2, Z.miss{p} must be a cell array; a plain logical array errors.
    rng(4, 'twister');
    K  = 5;
    sz = {10, 12 * ones(1, K), K};
    Z.model               = {'PAR2'};
    Z.loss_function       = {'Frobenius'};
    Z.loss_function_param = {[]};
    Z.modes               = {[1 2 3]};
    Z.size                = sz;
    Z.coupling.lin_coupled_modes   = [0 0 0];
    Z.coupling.coupling_type       = [];
    Z.coupling.coupl_trafo_matrices = cell(3,1);
    Z.constrained_modes   = [0 0 0];
    Z.constraints         = cell(3,1);
    Z.weights             = [1];
    for k = 1:K
        Z.object{1}{k} = randn(10, 12);
    end
    Z.miss{1} = true(10, 12);   % should be a cell, not a matrix

    init_options.lambdas_init = {ones(1, 2)};
    init_options.nvecs        = 0;
    init_options.distr        = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) rand(x,y)+0.1};
    init_options.normalize    = 1;
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

    testCase.verifyError(...
        @() cmtf.aoadmm.cmtf_AOADMM(Z, 'alg_options', fast_options(), ...
            'init', init_fac, 'init_options', init_options), ...
        'cmtf:missingData:PAR2maskNotCell', ...
        'PAR2 Z.miss that is not a cell should produce an error');
end
