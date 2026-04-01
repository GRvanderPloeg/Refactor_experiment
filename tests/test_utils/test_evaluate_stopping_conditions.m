function tests = test_evaluate_stopping_conditions
    tests = functiontests(localfunctions);
end

function opts = make_options(abs_tol, rel_tol)
    opts.AbsFuncTol  = abs_tol;
    opts.OuterRelTol = rel_tol;
end

function test_converged_by_abs(testCase)
% All four current values below AbsFuncTol → stop = true.
    opts = make_options(1e-4, 1e-6);
    stop = cmtf.utils.evaluate_stopping_conditions( ...
        1e-5, 1e-5, 1e-5, 1e-5, ...  % current
        1e-5, 1e-5, 1e-5, 1e-5, ...  % old (same → zero relative change)
        opts);
    verifyTrue(testCase, stop);
end

function test_converged_by_rel(testCase)
% Values above AbsFuncTol but relative change tiny → stop = true.
    opts = make_options(1e-10, 1e-4);
    f_old = 1.0;
    f_new = f_old * (1 - 1e-5); % relative change << OuterRelTol
    stop = cmtf.utils.evaluate_stopping_conditions( ...
        f_new, f_new, f_new, f_new, ...
        f_old, f_old, f_old, f_old, ...
        opts);
    verifyTrue(testCase, stop);
end

function test_not_converged(testCase)
% Large values with large relative change → stop = false.
    opts = make_options(1e-10, 1e-6);
    stop = cmtf.utils.evaluate_stopping_conditions( ...
        10, 10, 10, 10, ...
        20, 20, 20, 20, ...
        opts);
    verifyFalse(testCase, stop);
end

function test_partial_convergence(testCase)
% Three components converged, one (f_couplings) not → stop = false.
    opts = make_options(1e-10, 1e-4);
    f_old = 1.0;
    f_converged = f_old * (1 - 1e-5);
    f_not_converged = f_old * 0.5; % 50% relative change, above OuterRelTol
    stop = cmtf.utils.evaluate_stopping_conditions( ...
        f_converged, f_not_converged, f_converged, f_converged, ...
        f_old, f_old, f_old, f_old, ...
        opts);
    verifyFalse(testCase, stop);
end

function test_old_value_zero(testCase)
% When f_*_old == 0, the else branch computes absolute diff instead of ratio.
% With f_new also small (< AbsFuncTol), all components should converge.
    opts = make_options(1e-4, 1e-6);
    stop = cmtf.utils.evaluate_stopping_conditions( ...
        1e-5, 1e-5, 1e-5, 1e-5, ...  % current (below AbsFuncTol)
        0,    0,    0,    0,    ...  % old = 0 → triggers else branch
        opts);
    verifyTrue(testCase, stop);
end

function test_old_value_zero_not_converged(testCase)
% When f_*_old == 0 but f_new is large, the absolute diff exceeds OuterRelTol.
    opts = make_options(1e-10, 1e-6);
    stop = cmtf.utils.evaluate_stopping_conditions( ...
        1.0, 1.0, 1.0, 1.0, ...
        0,   0,   0,   0,   ...
        opts);
    verifyFalse(testCase, stop);
end
