% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: provide simple example of using the TrafficDensity class

% Instantiate object
td = TrafficDensity;

% Evaluate geographic area (fullest extent of data)
% Will use default plotting
td = td.run;