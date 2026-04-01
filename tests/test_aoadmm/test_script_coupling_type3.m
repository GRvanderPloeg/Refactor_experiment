% Test for coupling type 3 (C = H*Delta): binary subject selection
%
% Two CP tensors of order 3: sizes [20,30,40] and [15,25,35], both with
% 3 components. Modes 1 and 4 (the subject modes) are coupled via type 3.
%
% Delta is a global factor matrix of size K x R (K = 25 total unique
% subjects, R = 3 components). Binary H matrices select which subjects
% appear in each dataset:
%   H{1}: 20x25 - subjects 1-20 belong to dataset 1
%   H{4}: 15x25 - subjects 1-10 (shared) and 21-25 (dataset-2-only)
%
% The coupling constraint C{m} = H{m}*Delta enforces that factor rows for
% shared subjects (1-10) are identical across datasets, while subjects
% unique to one dataset (11-20 for dataset 1, 21-25 for dataset 2) are
% free to vary independently.
%
% Gaussian noise is added. The loss function is squared Frobenius norm.

function tests = test_script_coupling_type3
    tests = functiontests(localfunctions);
end

function test_workflow(testCase)

    % Control randomness
    rng(123, 'twister');

    % Specify synthetic data
    sz     = {20, 30, 40, 15, 25, 35}; % size of each mode
    P      = 2;                         % number of tensors
    lambdas_data = {[1 1 1], [1 1 1]};  % norms of components (R=3 each)
    modes  = {[1 2 3], [4 5 6]};        % modes per dataset
    noise  = 0.2;                       % noise level
    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    normalize_columns = 0;

    % Specify tensor model
    model{1} = 'CP';
    model{2} = 'CP';

    % Specify coupling: type 3 (C = H*Delta)
    % K = 25 total unique subjects:
    %   1-10:  shared (appear in both datasets)
    %   11-20: only in dataset 1
    %   21-25: only in dataset 2
    coupling.lin_coupled_modes    = [1 0 0 1 0 0]; % modes 1 and 4 are coupled
    coupling.coupling_type        = [3];            % coupling type 3
    coupling.coupl_trafo_matrices = cell(6, 1);

    % H{1}: 20x25, selects subjects 1-20 for dataset 1
    coupling.coupl_trafo_matrices{1} = zeros(20, 25);
    coupling.coupl_trafo_matrices{1}(1:20, 1:20) = eye(20);

    % H{4}: 15x25, selects subjects 1-10 (shared) and 21-25 (dataset-2-only)
    coupling.coupl_trafo_matrices{4} = zeros(15, 25);
    coupling.coupl_trafo_matrices{4}(1:10, 1:10) = eye(10);   % shared subjects
    coupling.coupl_trafo_matrices{4}(11:15, 21:25) = eye(5);  % dataset-2-only subjects

    % Loss functions
    loss_function{1}       = 'Frobenius';
    loss_function{2}       = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    %% Check model
    cmtf.utils.check_data_input(sz, modes, lambdas_data, coupling, loss_function, model);

    % Initialization options
    init_options.lambdas_init = {[1 1 1], [1 1 1]};
    init_options.nvecs        = 0;
    init_options.distr        = distr_data;
    init_options.normalize    = 1;

    % No constraints
    constrained_modes = [0 0 0 0 0 0];
    constraints       = cell(length(sz), 1);

    % Weights
    weights = [1/2 1/2];

    % Build Z struct
    Z.loss_function       = loss_function;
    Z.loss_function_param = loss_function_param;
    Z.model               = model;
    Z.modes               = modes;
    Z.size                = sz;
    Z.coupling            = coupling;
    Z.constrained_modes   = constrained_modes;
    Z.constraints         = constraints;
    Z.weights             = weights;

    % Create synthetic data
    [X, Atrue, ~, ~] = cmtf.utils.create_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, ...
        'lambdas', lambdas_data, 'noise', noise, ...
        'coupling', coupling, 'normalize_columns', normalize_columns, ...
        'distr_data', distr_data, 'loss_function', Z.loss_function);

    % Normalise objects
    normZ = cell(P, 1);
    for p = 1:P
        Z.object{p} = X{p};
        normZ{p}    = norm(Z.object{p});
        Z.object{p} = Z.object{p} / normZ{p};
    end

    % Create random initialization
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

    % Algorithm options
    options.Display                = 'no';
    options.DisplayIters           = 10;
    options.MaxOuterIters          = 20000; % This case converges slowly.
    options.MaxInnerIters          = 5;
    options.AbsFuncTol             = 1e-4;
    options.OuterRelTol            = 1e-8;
    options.innerRelPrTol_coupl    = 1e-3;
    options.innerRelPrTol_constr   = 1e-3;
    options.innerRelDualTol_coupl  = 1e-3;
    options.innerRelDualTol_constr = 1e-3;
    options.bsum                   = 0;
    options.eps_log                = 1e-10;

    % Run algorithm
    [Zhat, ~, ~, out] = cmtf.aoadmm.cmtf_AOADMM(Z, ...
        'alg_options', options, 'init', init_fac, 'init_options', init_options);

    % FIT
    Fit1 = 100 * (1 - norm(Z.object{1} - full(Zhat{1}))^2 / norm(Z.object{1})^2);
    Fit2 = 100 * (1 - norm(Z.object{2} - full(Zhat{2}))^2 / norm(Z.object{2})^2);

    % FMS
    true_ktensor{1} = ktensor(lambdas_data{1}' ./ normZ{1}, Atrue(modes{1}));
    FMS1 = score(Zhat{1}, true_ktensor{1});
    true_ktensor{2} = ktensor(lambdas_data{2}' ./ normZ{2}, Atrue(modes{2}));
    FMS2 = score(Zhat{2}, true_ktensor{2});

    % Assertions
    testCase.verifyEqual(out.f_tensors,       0.1690, "AbsTol", 1e-3);
    testCase.verifyEqual(Fit1,                89.0061, "AbsTol", 1e-3);
    testCase.verifyEqual(FMS1,                0.8767, "AbsTol", 1e-3);
    testCase.verifyEqual(Fit2,                77.1915, "AbsTol", 1e-3);
    testCase.verifyEqual(FMS2,                0.6959, "AbsTol", 1e-3);
    testCase.verifyEqual(out.OuterIterations,  17107);
end