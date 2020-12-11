% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%% Import plugins
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner

%% Run suite of tests
% https://www.mathworks.com/help/matlab/ref/runtests.html
% https://www.mathworks.com/help/matlab/matlab_prog/run-tests-for-various-workflows.html
testFolder = [getenv('TrafficDensityPath') filesep 'Testing'];
suite = testsuite(testFolder,'Name','UnitTest*');
runner = TestRunner.withNoPlugins;
results = runner.run(suite);

% Display status
if all([results.Passed])
    fprintf('Ran %i unit tests, all passed\n',numel(results));
else
    warning('At least one unit test failed');
end
%% Save
outFile = [getenv('TrafficDensityPath') filesep 'Testing' filesep datestr(now,'yyyymmddTHHMMSS') '_' version('-release') '_TestResults.mat'];
save(outFile,'suite','results');