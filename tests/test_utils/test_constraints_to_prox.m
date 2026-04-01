function tests = test_constraints_to_prox
    tests = functiontests(localfunctions);
end

% Helper: build minimal inputs for a single constrained mode.
% sz_val is the size used for mode 1.
function [constrained_modes, constraints, sz] = single_mode(constraint_spec, sz_val)
    constrained_modes = [1];
    constraints = {constraint_spec};
    sz = {sz_val};
end

function test_non_negativity(testCase)
    [cm, c, sz] = single_mode({'non-negativity'}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyEmpty(testCase, reg{1});
end

function test_box(testCase)
    [cm, c, sz] = single_mode({'box', 0, 1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_simplex_columnwise(testCase)
    [cm, c, sz] = single_mode({'simplex column-wise', 1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_simplex_rowwise(testCase)
    [cm, c, sz] = single_mode({'simplex row-wise', 1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_non_decreasing(testCase)
    [cm, c, sz] = single_mode({'non-decreasing'}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_non_increasing(testCase)
    [cm, c, sz] = single_mode({'non-increasing'}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_unimodality(testCase)
    [cm, c, sz] = single_mode({'unimodality', 1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_l1_ball(testCase)
    [cm, c, sz] = single_mode({'l1-ball', 1.0}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_l2_ball(testCase)
    [cm, c, sz] = single_mode({'l2-ball', 1.0}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_nonneg_l2_ball(testCase)
    [cm, c, sz] = single_mode({'non-negative l2-ball', 1.0}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_nonneg_l2_sphere(testCase)
    [cm, c, sz] = single_mode({'non-negative l2-sphere', 1.0}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_orthonormal(testCase)
    [cm, c, sz] = single_mode({'orthonormal'}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
end

function test_l1_regularization(testCase)
    [cm, c, sz] = single_mode({'l1 regularization', 0.1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_l0_regularization(testCase)
    [cm, c, sz] = single_mode({'l0 regularization', 0.1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_l2_regularization(testCase)
    [cm, c, sz] = single_mode({'l2 regularization', 0.1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_ridge(testCase)
    [cm, c, sz] = single_mode({'ridge', 0.1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_quadratic_regularization(testCase)
    n = 10;
    L = diag(2*ones(n,1)) + diag(-ones(n-1,1),1) + diag(-ones(n-1,1),-1);
    L(1,1) = 1; L(end,end) = 1;
    [cm, c, sz] = single_mode({'quadratic regularization', 0.1, L}, n);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_GL_smoothness(testCase)
    [cm, c, sz] = single_mode({'GL smoothness', 0.1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_TV_regularization(testCase)
    [cm, c, sz] = single_mode({'TV regularization', 0.1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_tPARAFAC2(testCase)
    [cm, c, sz] = single_mode({'tPARAFAC2', 0.1}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_custom_prox_only(testCase)
    my_prox = @(x, rho) max(x, 0);
    [cm, c, sz] = single_mode({'custom', my_prox}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyEmpty(testCase, reg{1});
end

function test_custom_with_reg_func(testCase)
    my_prox = @(x, rho) max(x, 0);
    my_reg  = @(x) 0;
    [cm, c, sz] = single_mode({'custom', my_prox, my_reg}, 10);
    [prox, reg] = cmtf.utils.constraints_to_prox(cm, c, sz);
    verifyTrue(testCase, isa(prox{1}, 'function_handle'));
    verifyTrue(testCase, isa(reg{1}, 'function_handle'));
end

function test_error_missing_constraint(testCase)
% A constrained mode with an empty constraint cell should error.
    constrained_modes = [1];
    constraints = {[]};
    sz = {10};
    verifyError(testCase, ...
        @() cmtf.utils.constraints_to_prox(constrained_modes, constraints, sz), ...
        '');
end

function test_unconstrained_mode_skipped(testCase)
% Mode with constrained_modes=0 should leave prox_operators{1} empty.
    constrained_modes = [0];
    constraints = {{'non-negativity'}};
    sz = {10};
    [prox, reg] = cmtf.utils.constraints_to_prox(constrained_modes, constraints, sz);
    verifyEmpty(testCase, prox{1});
    verifyEmpty(testCase, reg{1});
end
