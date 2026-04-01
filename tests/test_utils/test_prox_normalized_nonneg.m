function tests = test_prox_normalized_nonneg
    tests = functiontests(localfunctions);
end

function test_positive_column_unit_norm(testCase)
% A positive column [3;4] should normalise to [0.6;0.8].
    Y = cmtf.prox.prox_normalized_nonneg([3; 4]);
    verifyEqual(testCase, Y, [0.6; 0.8], 'AbsTol', 1e-14);
end

function test_negative_entries_zeroed(testCase)
% Negative entries must be zeroed before normalisation.
    Y = cmtf.prox.prox_normalized_nonneg([-1; 2; 3]);
    verifyEqual(testCase, Y(1), 0, 'AbsTol', 1e-14);
    verifyEqual(testCase, norm(Y, 2), 1, 'AbsTol', 1e-14);
end

function test_all_negative_column(testCase)
% All-negative column: entry with the largest value gets set to 1, rest 0.
% Input [-2;-1;-3]: max is -1 at row 2.
    Y = cmtf.prox.prox_normalized_nonneg([-2; -1; -3]);
    verifyEqual(testCase, Y(2), 1, 'AbsTol', 1e-14);
    verifyEqual(testCase, Y(1), 0, 'AbsTol', 1e-14);
    verifyEqual(testCase, Y(3), 0, 'AbsTol', 1e-14);
end

function test_all_zero_column(testCase)
% All-zero column: max coord (first, by MATLAB tie-breaking) gets set to 1.
    Y = cmtf.prox.prox_normalized_nonneg(zeros(4, 1));
    verifyEqual(testCase, sum(Y), 1, 'AbsTol', 1e-14);
    verifyEqual(testCase, norm(Y, 2), 1, 'AbsTol', 1e-14);
    verifyTrue(testCase, all(Y >= 0));
end
