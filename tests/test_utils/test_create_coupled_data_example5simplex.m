function tests = test_create_coupled_data_example5simplex
    tests = functiontests(localfunctions);
end

function test_CP_single(testCase)
    rng(0, 'twister');

    model{1} = 'CP';
    sz     = {60, 50, 40};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise  = 0.2;

    coupling.lin_coupled_modes    = [0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type0(testCase)
    rng(1, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20, 30, 40, 20, 70};
    modes  = {[1 2 3], [4 5]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.0;

    coupling.lin_coupled_modes    = [1 0 0 1 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(5,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_simplex_mode(testCase)
% Mode index 6 is uncoupled → columns of A{6} must sum to 1 (simplex constraint).
    rng(2, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    % Modes 1 and 4 are coupled (type 0); modes 5 and 6 are uncoupled mode-2/3 of tensor 2.
    sz     = {20, 30, 40, 20, 50, 15};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.0;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(6,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y), ...
                  @(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, A, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    % Each column of A{6} must sum to 1 (simplex).
    col_sums = sum(A{6}, 1);
    verifyEqual(testCase, col_sums, ones(1, length(lambdas_data{2})), 'AbsTol', 1e-12);
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
    d = 2;
    coupling.coupl_trafo_matrices{1} = rand(R, d);
    coupling.coupl_trafo_matrices{4} = rand(R, d);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type1(testCase)
    rng(4, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20, 30, 40, 40, 70, 80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 1;
    coupling.coupl_trafo_matrices = cell(6,1);
    coupling.coupl_trafo_matrices{1} = eye(20, 20);
    H = zeros(20, 40);
    for i = 1:20, H(i, 2*i-1) = 1; end
    coupling.coupl_trafo_matrices{4} = H;

    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type3(testCase)
    rng(5, 'twister');

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
    d = 5;
    coupling.coupl_trafo_matrices{1} = rand(sz{1}, d);
    coupling.coupl_trafo_matrices{4} = rand(sz{4}, d);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type4(testCase)
    rng(6, 'twister');

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
    coupling.coupl_trafo_matrices{4} = [eye(3,3); 0 0 0];

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_normalize_columns(testCase)
% normalize_columns=1: uncoupled factor matrix columns should have unit norm.
    rng(7, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20, 30, 40, 20, 70};
    modes  = {[1 2 3], [4 5]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(5,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) rand(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, A, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function, 'normalize_columns', 1);

    for r = 1:length(lambdas_data{1})
        verifyEqual(testCase, norm(A{2}(:,r)), 1, 'AbsTol', 1e-10);
        verifyEqual(testCase, norm(A{3}(:,r)), 1, 'AbsTol', 1e-10);
    end
end

function test_PAR2_single(testCase)
% Single PAR2 tensor: A{2} is a cell of shifted Bk matrices; X{1} is a cell.
    rng(8, 'twister');

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

    [X, A, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    verifyTrue(testCase, iscell(A{2}));
    verifyEqual(testCase, length(A{2}), 10);
    verifyTrue(testCase, iscell(X{1}));
end

function test_KL_loss(testCase)
% Single CP tensor with KL loss: output must contain non-negative integers.
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
    loss_function = {'KL'};

    [X, ~, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    verifyTrue(testCase, all(double(X{1}(:)) >= 0));
    verifyTrue(testCase, all(double(X{1}(:)) == floor(double(X{1}(:)))));
end

function test_IS_loss(testCase)
% Single CP tensor with IS (gamma) loss: output must be positive.
    rng(10, 'twister');

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

    [X, ~, ~, ~] = cmtf.utils.create_coupled_data_example5simplex( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function, 'loss_function_param', loss_function_param);

    verifyTrue(testCase, all(double(X{1}(:)) > 0));
end
