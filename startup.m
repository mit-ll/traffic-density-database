% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Script will add necessary paths

[~,currDir] = fileparts(pwd);
if ~strcmp(currDir,'TrafficDensityDatabase')
    error('Startup must be executed from TrafficDensityDatabase directory');
end

disp('Running TrafficDensityDatabase startup script...');

disp('Adding paths...')
addpath('.');
addpath(genpath('Data'));
addpath(genpath('Utilities'));
addpath(genpath('Examples'));
addpath(genpath('Testing'));

setenv('TrafficDensityPath',pwd);

% Try to build necessary files if they do not exist
if ~exist('accumarraymax_mex','file')
    disp('Building accumarraymax_mex...')
    cd('Utilities');
    buildaccumarray;
    cd('..');    
end
    
% Verify correct version of Matlab
disp('Checking Matlab version...')
if verLessThan('matlab','9.5')
    warning('Software only tested on Matlab R2018b and newer')
end

% Ensure that have necessary toolboxes to run software
disp('Checking licenses...')
reqlicenses = {'map_toolbox','matlab','statistics_toolbox'};
insufLicense = false;
for rr = 1:length(reqlicenses)
    if ~license('test',reqlicenses{rr})
        warning('Must have %s toolbox installed to use all features of the traffic density tool',reqlicenses{rr})
        insufLicense = true;
    end
end
if insufLicense
    error('One or more required toolboxes is not installed (see above warning for specific toolbox(es))');
end

% Everything is done
disp('Done!');