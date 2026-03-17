function tests = test_create_coupled_data_smoothBks
    tests = functiontests(localfunctions);
end

function [X, A] = generate_data(seed)
    rng(seed, 'twister');

    model{1} = 'CP';
    model{2} = 'PAR2';
    sz     = {30, 20, 50, 30, 200*ones(1,10), 10};
    modes  = {[1 2 3], [4 5 6]};
    lambdas_data = {[1 1 1], [1 1 1]};
    noise  = 0.5;

    coupling.lin_coupled_modes    = [1 0 0 1 0 0];
    coupling.coupling_type        = 0;
    coupling.coupl_trafo_matrices = cell(6,1);

    normalize_columns = 0;
    distr_data = {@(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) randn(x,y), ...
                  @(x,y) randn(x,y), @(x,y) randn(x,y), @(x,y) rand(x,y)+0.1};

    loss_function = {'Frobenius','Frobenius'};

    [X, A, ~, ~] = cmtf.utils.create_coupled_data_smoothBks( ...
        'model', model, 'size', sz, 'modes', modes, 'lambdas', lambdas_data, ...
        'noise', noise, 'coupling', coupling, 'normalize_columns', normalize_columns, ...
        'distr_data', distr_data, 'loss_function', loss_function);
end

function test_reproducibility(testCase)
% Two calls with the same seed must produce bit-identical results.
% This guards against LAPACK-internal non-determinism (e.g. DGESDD random
% vectors) or platform-dependent QR/SVD sign conventions.
    [X1, A1] = generate_data(42);
    [X2, A2] = generate_data(42);

    % CP tensor (mode 1): must be identical
    verifyEqual(testCase, double(X1{1}), double(X2{1}));

    % PAR2 slices (mode 2): must be identical
    for k = 1:length(X1{2})
        verifyEqual(testCase, X1{2}{k}, X2{2}{k});
    end

    % Smooth Bk factor matrices must be identical
    for k = 1:length(A1{5})
        verifyEqual(testCase, A1{5}{k}, A2{5}{k});
    end
end

function test_smoothBks_lie_in_polynomial_subspace(testCase)
% Each Bk must lie in the column space of the polynomial basis
% M = [1, x, x^2, x^3] evaluated on linspace(-1,1,200).
% Equivalently, projecting Bk onto the orthogonal complement of that
% subspace should give (near) zero.
    [~, A] = generate_data(42);

    szBk = 200;
    x = linspace(-1, 1, szBk);
    M = [ones(szBk,1), x', x.^2', x.^3'];
    [Morth, ~] = qr(M, 0);
    P_perp = eye(szBk) - Morth*Morth'; % projector onto orthogonal complement

    for k = 1:length(A{5})
        residual = norm(P_perp * A{5}{k}, 'fro');
        verifyLessThan(testCase, residual, 1e-10);
    end
end
