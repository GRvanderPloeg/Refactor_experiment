function tests = test_project_ortho
    tests = functiontests(localfunctions);
end

function test_tall_matrix_orthonormal(testCase)
% Tall matrix (m>n): columns of output must be orthonormal.
    rng(0, 'twister');
    Z = cmtf.project.project_ortho(rand(10, 3));
    verifyEqual(testCase, Z'*Z, eye(3), 'AbsTol', 1e-12);
end

function test_square_matrix_orthonormal(testCase)
% Square matrix: output must be orthonormal (Z'*Z = I).
    rng(1, 'twister');
    Z = cmtf.project.project_ortho(rand(5, 5));
    verifyEqual(testCase, Z'*Z, eye(5), 'AbsTol', 1e-12);
end

function test_output_size(testCase)
% Output size must equal input size.
    rng(2, 'twister');
    X = rand(8, 4);
    Z = cmtf.project.project_ortho(X);
    verifyEqual(testCase, size(Z), size(X));
end

function test_already_orthonormal(testCase)
% An already-orthonormal input should remain orthonormal.
    rng(3, 'twister');
    Q = orth(rand(6, 3));
    Z = cmtf.project.project_ortho(Q);
    verifyEqual(testCase, Z'*Z, eye(3), 'AbsTol', 1e-12);
end
