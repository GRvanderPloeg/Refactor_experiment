function tests = test_create_coupled_data
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Runs ONCE for this file, before all tests here
    setup;       % your project's startup, e.g. add paths
    rng(0,'twister');
end

function setup(testCase)
    % Runs before EACH test function in this file
    rng(0,'twister');  % if you want deterministic randomness per test
end

function test_CP_CP(testCase)

    model{1} = 'CP';
    model{2} = 'CP';
    sz     = {20,30,40, 20,30,20};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data= {[1 1 1], [1 1 1]};
    noise = 0;

    coupling.lin_coupled_modes = [1 0 0 1 0 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(6,1);
    
    normalize_columns = 0;
    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    [~, ~, ~,~] = cmtf.utils.create_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',loss_function); %create data
end

function test_CP_PAR2(testCase)

    model{1} = 'CP';
    model{2} = 'PAR2';
    sz     = {20,30,40, 20,30*ones(1,20),20};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data= {[1 1 1], [1 1 1]};
    noise = 0;

    coupling.lin_coupled_modes = [1 0 0 1 0 0];
    coupling.coupling_type = 0;
    coupling.coupl_trafo_matrices = cell(6,1);
    
    normalize_columns = 0;
    distr_data = {@(x,y) rand(x,y), @(x,y) randn(x,y),@(x,y) randn(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y),@(x,y) rand(x,y)+0.1};

    loss_function{1} = 'Frobenius';
    loss_function{2} = 'Frobenius';

    [~, ~, ~,~] = cmtf.utils.create_coupled_data('model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, 'noise', noise,'coupling',coupling,'normalize_columns',normalize_columns,'distr_data',distr_data,'loss_function',loss_function); %create data
end