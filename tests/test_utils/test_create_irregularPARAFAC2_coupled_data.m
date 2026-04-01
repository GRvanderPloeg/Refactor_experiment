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
