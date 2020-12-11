% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: provide simple example of using the TrafficDensity class
% with lowh cell (up to 5000 ft AGL and finer horizontal discretization - 2.5 NM)
td = TrafficDensity('joblowh.dat');

% Specify low-altitude tables
td.filenames.cell = 'celllowh.mat';
td.filenames.cellCoverage = 'cellcoveragelowh.mat';
td.filenames.cellAirspace = 'cellAirspacelowh.mat';

td = td.LoadData;

td = td.run;
     