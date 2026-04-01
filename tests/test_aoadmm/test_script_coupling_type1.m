% Test for coupling type 1 (H*C = Delta): binary subject selection
%
% Two CP tensors of order 3: sizes [20,30,40] and [15,25,35], both with
% 3 components. Modes 1 and 4 (the subject modes) are coupled via type 1.
%
% Binary H matrices specify which subjects are shared between datasets:
%   H{1}: 10x20 - selects 10 out of 20 subjects from tensor 1
%                  (subjects 2,4,6,8,10,12,14,16,18,20)
%   H{4}: 10x15 - selects 10 out of 15 subjects from tensor 2
%                  (subjects 1,2,3,4,5,6,7,8,9,10)
%
% The coupling constraint H{1}*C{1} = Delta = H{4}*C{4} enforces that the
% selected (shared) subjects have matching factor values across datasets.
% Subjects not selected by H are uncoupled and free to vary independently.
%
% Gaussian noise is added. The loss function is squared Frobenius norm.

function tests = test_script_coupling_type1
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

    % Specify coupling: type 1 (H*C = Delta)
    coupling.lin_coupled_modes    = [1 0 0 1 0 0]; % modes 1 and 4 are coupled
    coupling.coupling_type        = [1];            % coupling type 1
    coupling.coupl_trafo_matrices = cell(6, 1);

    % Binary H matrices for subject selection
    % H{1}: 10x20, selects every second subject from tensor 1
    coupling.coupl_trafo_matrices{1} = zeros(10, 20);
    for i = 1:10
        coupling.coupl_trafo_matrices{1}(i, 2*i) = 1;
    end
    % H{4}: 10x15, selects the first 10 subjects from tensor 2
    coupling.coupl_trafo_matrices{4} = zeros(10, 15);
    for i = 1:10
        coupling.coupl_trafo_matrices{4}(i, i) = 1;
    end

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
    options.MaxOuterIters          = 4000;
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
    testCase.verifyEqual(out.f_tensors,       0.0380, "AbsTol", 1e-3);
    testCase.verifyEqual(Fit1,                96.1944, "AbsTol", 1e-3);
    testCase.verifyEqual(FMS1,                0.9992, "AbsTol", 1e-3);
    testCase.verifyEqual(Fit2,                96.2115, "AbsTol", 1e-3);
    testCase.verifyEqual(FMS2,                0.9985, "AbsTol", 1e-3);
    testCase.verifyEqual(out.OuterIterations,  16);
end
