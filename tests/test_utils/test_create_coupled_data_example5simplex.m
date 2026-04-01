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
