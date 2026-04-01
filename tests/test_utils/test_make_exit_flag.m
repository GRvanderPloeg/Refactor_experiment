function tests = test_make_exit_flag
    tests = functiontests(localfunctions);
end

function opts = make_options(max_iters, abs_tol)
    opts.MaxOuterIters = max_iters;
    opts.AbsFuncTol    = abs_tol;
end

function test_max_iterations(testCase)
% iter > MaxOuterIters → 'maxIterations'.
    opts = make_options(100, 1e-6);
    flag = cmtf.utils.make_exit_flag(101, 0.5, 0.5, 0.5, 0.5, opts, false);
    verifyEqual(testCase, flag, 'maxIterations');
end

function test_illconditioned(testCase)
% illconditioned=true → "illconditioned lin system".
    opts = make_options(100, 1e-6);
    flag = cmtf.utils.make_exit_flag(50, 0.5, 0.5, 0.5, 0.5, opts, true);
    verifyEqual(testCase, flag, "illconditioned lin system");
end

function test_all_absfunctol(testCase)
% All f_* below AbsFuncTol → all struct fields "AbsFuncTol".
    opts = make_options(100, 1e-3);
    flag = cmtf.utils.make_exit_flag(50, 1e-4, 1e-4, 1e-4, 1e-4, opts, false);
    verifyEqual(testCase, flag.f_tensors,        "AbsFuncTol");
    verifyEqual(testCase, flag.f_couplings,      "AbsFuncTol");
    verifyEqual(testCase, flag.f_constraints,    "AbsFuncTol");
    verifyEqual(testCase, flag.f_PAR2_couplings, "AbsFuncTol");
end

function test_all_relfunctol(testCase)
% All f_* above AbsFuncTol → all struct fields "RelFuncTol".
    opts = make_options(100, 1e-10);
    flag = cmtf.utils.make_exit_flag(50, 1.0, 1.0, 1.0, 1.0, opts, false);
    verifyEqual(testCase, flag.f_tensors,        "RelFuncTol");
    verifyEqual(testCase, flag.f_couplings,      "RelFuncTol");
    verifyEqual(testCase, flag.f_constraints,    "RelFuncTol");
    verifyEqual(testCase, flag.f_PAR2_couplings, "RelFuncTol");
end

function test_mixed(testCase)
% f_tensors below AbsFuncTol, others above → mixed fields.
    opts = make_options(100, 1e-3);
    flag = cmtf.utils.make_exit_flag(50, 1e-4, 1.0, 1.0, 1.0, opts, false);
    verifyEqual(testCase, flag.f_tensors,        "AbsFuncTol");
    verifyEqual(testCase, flag.f_couplings,      "RelFuncTol");
    verifyEqual(testCase, flag.f_constraints,    "RelFuncTol");
    verifyEqual(testCase, flag.f_PAR2_couplings, "RelFuncTol");
end
