function tests = test_create_irregularPARAFAC2_coupled_data
    tests = functiontests(localfunctions);
end

function test_irregular_PAR2_single(testCase)
% Single irregular PAR2: verify A{2} is a cell of matrices.
    rng(0, 'twister');

    model{1} = 'PAR2';
    sz     = {40, [61:1:80], 20};  % 20 slices with varying sizes 61..80
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise  = 0.2;

    coupling.lin_coupled_modes    = [0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) rand(x,y)+0.1};
    loss_function = {'Frobenius'};

    [~, A, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    verifyTrue(testCase, iscell(A{2}));
    verifyEqual(testCase, length(A{2}), length(sz{2}));
end

function test_orthonormality(testCase)
% Each Bk must satisfy Bk'*Bk ≈ I (orthonormal columns).
    rng(1, 'twister');

    R = 3;
    model{1} = 'PAR2';
    sz     = {40, [50:1:59], 10};  % 10 slices with sizes 50..59
    modes  = {[1 2 3]};
    lambdas_data = {ones(1,R)};
    noise  = 0.0;

    coupling.lin_coupled_modes    = [0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) rand(x,y)+0.1};
    loss_function = {'Frobenius'};

    [~, A, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    for k = 1:length(sz{2})
        Bk = A{2}{k};
        verifyEqual(testCase, Bk' * Bk, eye(R), 'AbsTol', 1e-10);
    end
end

function test_matrix_PAR2_coupled(testCase)
% Matrix + irregular PAR2 coupled via coupling type 0. Exercises Delta generation.
    rng(2, 'twister');

    model{1} = 'CP';
    model{2} = 'PAR2';
    sz     = {40, 60, 40, [61:1:70], 10};  % 10 irregular slices
    modes  = {[1 2], [3 4 5]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.3;

    coupling.lin_coupled_modes    = [1 0 1 0 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(5,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) rand(x,y), ...
                  @(x,y) randn(x,y), @(x,y) rand(x,y)+0.1};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type1(testCase)
% Two CP tensors coupled via type 1 (linear transformation).
    rng(3, 'twister');

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

    [~, ~, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type2(testCase)
% Two CP tensors coupled via type 2 (component-space coupling).
    rng(4, 'twister');

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

    [~, ~, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type3(testCase)
% Two CP tensors coupled via type 3 (latent-space coupling).
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

    [~, ~, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_CP_CP_coupling_type4(testCase)
% Two CP tensors coupled via type 4 (partial coupling).
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

    [~, ~, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end

function test_normalize_columns(testCase)
% Coupling type 0 with normalize_columns=1: uncoupled factor columns have unit norm.
    rng(7, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20, 30, 40, 20, 70, 80};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.1;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(6,1);

    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y)};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, A, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function, 'normalize_columns', 1);

    for r = 1:length(lambdas_data{1})
        verifyEqual(testCase, norm(A{2}(:,r)), 1, 'AbsTol', 1e-10);
        verifyEqual(testCase, norm(A{3}(:,r)), 1, 'AbsTol', 1e-10);
    end
end

function test_KL_and_IS_loss(testCase)
% Two CP tensors: one with KL loss (Poisson), one with IS loss (gamma).
    rng(8, 'twister');

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20, 15, 10, 20, 12, 8};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.0;

    coupling.lin_coupled_modes    = [0 0 0 0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(6,1);

    distr_data = {@(x,y) rand(x,y)+0.5, @(x,y) rand(x,y)+0.5, @(x,y) rand(x,y)+0.5, ...
                  @(x,y) rand(x,y)+0.5, @(x,y) rand(x,y)+0.5, @(x,y) rand(x,y)+0.5};
    loss_function = {'KL', 'IS'};
    loss_function_param = {[], 2};

    [X, ~, ~, ~] = cmtf.utils.create_irregularPARAFAC2_coupled_data( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function, 'loss_function_param', loss_function_param);

    % KL: non-negative integers
    verifyTrue(testCase, all(double(X{1}(:)) >= 0));
    verifyTrue(testCase, all(double(X{1}(:)) == floor(double(X{1}(:)))));
    % IS: positive values
    verifyTrue(testCase, all(double(X{2}(:)) > 0));
end
