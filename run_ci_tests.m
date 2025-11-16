function run_ci_tests
%RUN_CI_TESTS Entry point for GitHub Actions CI.
%
% This function is called by the GitHub Actions workflow to:
%   1) Run startup.m to set paths
%   2) Run all tests in tests/ (including subfolders)
%   3) Fail CI if any tests fail

    % 1. Set up paths etc.
    if exist('setup', 'file')
        setup;
    else
        warning('startup.m not found on path.');
    end

    % 2. Run tests
    results = runtests('tests', 'IncludeSubfolders', true);

    % Optional: print a nice summary
    disp(table(results));

    % 3. Fail the CI run if any test fails
    assertSuccess(results);
end
