% Test for coupling type 4 (C = Delta*H): shared and distinct components
%
% One CP tensor [50,30,40] with R1=4 components and one CP matrix [50,70]
% with R2=3 components. Modes 1 and 4 (both size 50) are coupled via
% type 4.
%
% Delta is I x r (50 x 5), a global factor matrix with r=5 latent
% components. Binary H matrices assign each dataset component to exactly
% one Delta column:
%   H{1}: 5x4 - maps Delta columns 1-4 to tensor components 1-4
%   H{4}: 5x3 - maps Delta columns 1,2,5 to matrix components 1,2,3
%
% This gives:
%   Components 1,2: shared (both datasets use Delta columns 1,2)
%   Components 3,4 of tensor: distinct (only tensor uses Delta columns 3,4)
%   Component 3 of matrix:   distinct (only matrix uses Delta column 5)
%
% Gaussian noise is added. The loss function is squared Frobenius norm.

function tests = test_script_coupling_type4
    tests = functiontests(localfunctions);
end

function test_workflow(testCase)

    % Control randomness
    rng(123, 'twister');

    % Specify synthetic data
    sz           = {50, 30, 40, 50, 70};
    P            = 2;
    lambdas_data = {[1 1 1 1], [1 1 1]};   % R1=4, R2=3
    modes        = {[1 2 3], [4 5]};
    noise        = 0.2;
    distr_data   = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                    @(x,y) rand(x,y), @(x,y) rand(x,y)};
    normalize_columns = 0;

    % Specify tensor model
    model{1} = 'CP';
    model{2} = 'CP';

    % Specify coupling (type 4: C = Delta*H)
    % r = 5 total latent components:
    %   1,2:   shared between tensor and matrix
    %   3,4:   distinct to tensor (dataset 1)
    %   5:     distinct to matrix (dataset 2)
    coupling.lin_coupled_modes    = [1 0 0 1 0];
    coupling.coupling_type        = [4];
    coupling.coupl_trafo_matrices = cell(5, 1);

    % H{1}: 5x4, maps Delta columns 1-4 to tensor components 1-4
    coupling.coupl_trafo_matrices{1} = [eye(4); zeros(1, 4)];  % 5x4

    % H{4}: 5x3, maps Delta columns 1,2,5 to matrix components 1,2,3
    coupling.coupl_trafo_matrices{4} = zeros(5, 3);
    coupling.coupl_trafo_matrices{4}(1, 1) = 1;  % shared comp 1
    coupling.coupl_trafo_matrices{4}(2, 2) = 1;  % shared comp 2
    coupling.coupl_trafo_matrices{4}(5, 3) = 1;  % distinct comp 3

    % Loss functions
    loss_function{1}       = 'Frobenius';
    loss_function{2}       = 'Frobenius';
    loss_function_param{1} = [];
    loss_function_param{2} = [];

    % Check model validity
    cmtf.utils.check_data_input(sz, modes, lambdas_data, coupling, loss_function, model);

    % Initialization options
    init_options.lambdas_init = lambdas_data;
    init_options.nvecs        = 0;
    init_options.distr        = distr_data;
    init_options.normalize    = 1;

    % No constraints
    constrained_modes = [0 0 0 0 0];
    constraints       = cell(length(sz), 1);

    % Weights
    weights = [0.5 0.5];

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

    % Initialise factors
    init_fac = cmtf.utils.init_coupled_AOADMM_CMTF(Z, 'init_options', init_options);

    % Algorithm options
    options.Display               = 'iter';
    options.DisplayIters          = 10;
    options.MaxOuterIters         = 4000;
    options.MaxInnerIters         = 5;
    options.AbsFuncTol            = 1e-4;
    options.OuterRelTol           = 1e-8;
    options.innerRelPrTol_coupl   = 1e-3;
    options.innerRelPrTol_constr  = 1e-3;
    options.innerRelDualTol_coupl = 1e-3;
    options.innerRelDualTol_constr= 1e-3;
    options.bsum                  = 0;
    options.eps_log               = 1e-10;

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
    testCase.verifyEqual(out.f_tensors,       0.044979, "AbsTol", 1e-3);
    testCase.verifyEqual(Fit1,                96.1826, "AbsTol", 1e-3);
    testCase.verifyEqual(FMS1,                0.9975, "AbsTol", 1e-3);
    testCase.verifyEqual(Fit2,                94.8216, "AbsTol", 1e-3);
    testCase.verifyEqual(FMS2,                0.6154, "AbsTol", 1e-3);
    testCase.verifyEqual(out.OuterIterations,  27);
end