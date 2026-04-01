function tests = test_create_CP_data_example10piecewiseconstant
    tests = functiontests(localfunctions);
end

function test_CP_single(testCase)
% Single CP tensor: verify piecewise constant structure in mode-1 factor.
    rng(0, 'twister');

    model{1} = 'CP';
    sz     = {50, 30, 40};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise  = 0.1;

    coupling.lin_coupled_modes  = [0 0 0];
    coupling.coupling_type      = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius'};

    [~, A, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    % Each column of A{1} must have at most 5 distinct values (piecewise constant).
    for r = 1:length(lambdas_data{1})
        verifyLessThanOrEqual(testCase, numel(unique(A{1}(:,r))), 5);
    end
end

function test_CP_CP_coupling_type0(testCase)
    rng(1, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {50, 30, 40, 50, 70, 80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(6,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type1(testCase)
    rng(2, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {50, 30, 40, 100, 70, 80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 1;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = eye(50, 50);
    % Downsampling matrix: map 100-dim to 50-dim by taking every other entry
    H = zeros(50, 100);
    for i = 1:50, H(i, 2*i-1) = 1; end
    coupling.coupl_trafo_matrices{4} = H;

    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type2(testCase)
    rng(3, 'twister');

    R = 3;
    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20, 30, 40, 20, 70, 80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,R), ones(1,R)};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 2;
    coupling.coupl_trafo_matrices = cell(6,1);
    d = 2; % d <= R
    coupling.coupl_trafo_matrices{1} = rand(R, d);
    coupling.coupl_trafo_matrices{4} = rand(R, d);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type3(testCase)
    rng(4, 'twister');

    R = 3;
    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20, 30, 40, 25, 70, 80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {ones(1,R), ones(1,R)};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 3;
    coupling.coupl_trafo_matrices = cell(6,1);
    d = 5; % shared latent dimension
    coupling.coupl_trafo_matrices{1} = rand(sz{1}, d);
    coupling.coupl_trafo_matrices{4} = rand(sz{4}, d);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type4(testCase)
    rng(5, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {50, 30, 40, 50, 70};
    modes  = {[1 2 3], [4 5]};
    lambdas_data = {[1 1 1 1], [1 1 1]};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0];
    coupling.coupling_type        = 4;
    coupling.coupl_trafo_matrices = cell(5,1);
    coupling.coupl_trafo_matrices{1} = eye(4, 4);
    coupling.coupl_trafo_matrices{4} = [eye(3,3); 0 0 0]; % 4x3

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_normalize_columns(testCase)
% Coupling type 0 with normalize_columns=1: uncoupled factor columns have unit norm.
    rng(6, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {50, 30, 40, 50, 70, 80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(6,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, A, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function, 'normalize_columns', 1);

    for r = 1:length(lambdas_data{1})
        verifyEqual(testCase, norm(A{2}(:,r)), 1, 'AbsTol', 1e-10);
        verifyEqual(testCase, norm(A{3}(:,r)), 1, 'AbsTol', 1e-10);
    end
end

function test_PAR2_single(testCase)
% Single PAR2 tensor: A{2} is a cell of Bk matrices (SHIFT PARAFAC); X{1} is a cell.
    rng(7, 'twister');

    model{1} = 'PAR2';
    sz     = {40, 120*ones(1,10), 10};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise  = 0.2;

    coupling.lin_coupled_modes    = [0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y)+0.1};
    loss_function = {'Frobenius'};

    [X, A, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    verifyTrue(testCase, iscell(A{2}));
    verifyEqual(testCase, length(A{2}), 10);
    verifyTrue(testCase, iscell(X{1}));
end

function test_KL_loss(testCase)
% Single CP tensor with KL (Poisson) loss.
% Note: mode-1 factor is always piecewise-constant in [-1,1] so the tensor
% can have negative entries; poissrnd returns NaN for those.  The purpose of
% this test is to exercise the KL code path, so we verify the output size only.
    rng(8, 'twister');

    model{1} = 'CP';
    sz     = {20, 15, 10};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise  = 0.0;

    coupling.lin_coupled_modes    = [0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) rand(x,y)+0.5, @(x,y) rand(x,y)+0.5, @(x,y) rand(x,y)+0.5};
    loss_function = {'KL'};

    [X, ~, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    verifyEqual(testCase, length(X), 1);
    verifyTrue(testCase, isa(X{1}, 'tensor'));
end

function test_IS_loss(testCase)
% Single CP tensor with IS (gamma) loss.
% Note: mode-1 factor is always piecewise-constant in [-1,1] so the tensor
% can have negative entries.  The purpose of this test is to exercise the IS
% code path, so we verify the output size only.
    rng(9, 'twister');

    model{1} = 'CP';
    sz     = {20, 15, 10};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise  = 0.0;

    coupling.lin_coupled_modes    = [0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) rand(x,y)+0.5, @(x,y) rand(x,y)+0.5, @(x,y) rand(x,y)+0.5};
    loss_function = {'IS'};
    loss_function_param = {2};

    [X, ~, ~, ~] = cmtf.utils.create_CP_data_example10piecewiseconstant( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function, 'loss_function_param', loss_function_param);

    verifyEqual(testCase, length(X), 1);
    verifyTrue(testCase, isa(X{1}, 'tensor'));
end
