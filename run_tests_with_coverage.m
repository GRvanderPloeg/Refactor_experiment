function run_tests_with_coverage()
%RUN_TESTS_WITH_COVERAGE  Run all tests and generate a Cobertura coverage report.
%   Produces coverage.xml in the project root covering the +cmtf package.

    % Ensure dependencies are on the path
    setup();

    % Discover all tests
    suite = matlab.unittest.TestSuite.fromFolder('tests', 'IncludingSubfolders', true);

    % Set up runner with text output and coverage
    runner = matlab.unittest.TestRunner.withTextOutput;
    rootDir = fileparts(mfilename('fullpath'));
    cmtfDir = fullfile(rootDir, '+cmtf');
    plugin = matlab.unittest.plugins.CodeCoveragePlugin.forFolder(cmtfDir, ...
        'IncludingSubfolders', true, ...
        'Producing', matlab.unittest.plugins.codecoverage.CoberturaFormat('coverage.xml'));
    runner.addPlugin(plugin);

    % Run
    results = runner.run(suite);

    % Print summary
    disp(results.table);

    % Fail loudly if any tests failed
    assertSuccess(results);
end