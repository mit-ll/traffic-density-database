% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: provide simple example of using the TrafficDensity class, with
% the user specifying a geographic region

% Instantiate object
td = TrafficDensity;

% Evaluate geographic area (based on user specified limits) Data is gridded
% to 1/12 degree by default (1/24 degree below 5000 ft). 1/48 degree is
% used to test Matlab's binning functionality.
td.area.LatitudeLimit = [35+1/48,40];
td.area.LongitudeLimit = [-90,-85-1/48]; %Note: Negative numbers correspond to West Longitude (i.e., -90 --> 90 W)
td = td.run;
