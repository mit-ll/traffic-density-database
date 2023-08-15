% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: provide simple example of using the TrafficDensity class, with
% the user specifying a track

% Instantiate object
td = TrafficDensity;

% Load example track
t = readtable('exampleTrack.csv');

computeinds = 1:15:length(t.time); % Reduce number of points to speed up processing - recommended time between updates is 5-30 s

% Set properties in object
td.track.Time_s = t.time(computeinds);
td.track.Latitude_deg = t.lat(computeinds);
td.track.Longitude_deg = t.lon(computeinds);
td.track.Altitude_MSL_ft = t.altitude(computeinds);
td.track.Speed_kts = t.tas(computeinds);

% Compute standard deviation
td.computestd = true;

% Evaluate track collision risk
td = td.run;
     