function tests = test_create_coupled_data_unimodalBks
    tests = functiontests(localfunctions);
end

function test_PAR2_single(testCase)
% Single regular PAR2 tensor with unimodal Bk matrices.
    rng(0, 'twister');

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

    [X, A, ~, ~] = cmtf.utils.create_coupled_data_unimodalBks( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    verifyEqual(testCase, length(X), 1);
    verifyEqual(testCase, length(A{2}), 10); % cell of Bk matrices
end

function test_unimodal_shape(testCase)
% Each column of A{2}{1} should be unimodal: diff changes sign at most once
% (from positive to negative).
    rng(1, 'twister');

    model{1} = 'PAR2';
    sz     = {40, 120*ones(1,5), 5};
    modes  = {[1 2 3]};
    lambdas_data = {[1 1 1]};
    noise  = 0.0;

    coupling.lin_coupled_modes    = [0 0 0];
    coupling.coupling_type        = [];
    coupling.coupl_trafo_matrices = cell(3,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) rand(x,y), @(x,y) rand(x,y)+0.1};
    loss_function = {'Frobenius'};

    [~, A, ~, ~] = cmtf.utils.create_coupled_data_unimodalBks( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);

    R = length(lambdas_data{1});
    for r = 1:R
        col = A{2}{1}(:,r);
        d   = diff(col);
        % Sign changes from + to - count as peaks; each column should have at most 1 peak.
        sign_changes = sum(diff(sign(d(d~=0))) < 0);
        verifyLessThanOrEqual(testCase, sign_changes, 1);
    end
end

function test_matrix_PAR2_coupled(testCase)
% Matrix + regular PAR2 coupled via coupling type 0.
    rng(2, 'twister');

    model{1} = 'CP';
    model{2} = 'PAR2';
    sz     = {40, 60, 40, 120*ones(1,20), 20};
    modes  = {[1 2], [3 4 5]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.3;

    coupling.lin_coupled_modes    = [1 0 1 0 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(5,1);

    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y), @(x,y) rand(x,y), ...
                  @(x,y) rand(x,y), @(x,y) rand(x,y)+0.1};
    loss_function = {'Frobenius', 'Frobenius'};

    [~, ~, ~, ~] = cmtf.utils.create_coupled_data_unimodalBks( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'distr_data', distr_data, ...
        'loss_function', loss_function);
end
