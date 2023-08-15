% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: provide simple example of using the TrafficDensity class, with
% the user specifying temporal characteristics

% Instantiate object
td = TrafficDensity;

% Specify certain temporal characteristics
td.timeofday = [2,3];       % Approximately night (UTC)
td.monthofyear = [6,7,8];   % Summer months
td.dayofweek = [2:6];       % Weekdays

% Specify only certain aircraft types of interest (1200-code)
td.ACcategory = 1; % 1 - 1200-code, 0 - discrete-code

% Load track
t = readtable('exampleTrack.csv');
computeinds = 1:15:length(t.time); % Reduce number of points to speed up processing - recommended time between updates is 5-30 s
td.track.Time_s = t.time(computeinds);
td.track.Latitude_deg = t.lat(computeinds);
td.track.Longitude_deg = t.lon(computeinds);
td.track.Altitude_MSL_ft = t.altitude(computeinds);
td.track.Speed_kts = t.tas(computeinds);

% Execute processing
td = td.run;
     