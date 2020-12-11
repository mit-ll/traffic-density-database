% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: provide example of using the TrafficDensity class, with plotting
% options, also add noncooperative intruder aircraft

% Instantiate object
td = TrafficDensity;

% Evaluate area (based on user specified limits)
td.area.LatitudeLimit = [37,41];
td.area.LongitudeLimit = [-(109+3/60),-(102+3/60)];

td.plotresults = false; % Save plotting for below

td.height = [0,1,2,3]; % Process only the first 4 altitude bins

td = td.run; % Execute processing

td.processNoncoop = true; % Add plotting of noncooperatives (it is only the plotting function that considers noncooperatives)

% Plot coverage information
td.plot('plottype','coverage');

% Plot airspace class
td.plot('plottype','airspace');

% Plot density with all options (see plot method documentation for other parameters)
% The outdata and summarize outputs provide user access to the estimated parameters
[outdata,summarize] = td.plot('plottype','density',...
    'plotscale','linear',...
    'plotvalue','avgupperbound',...
    'plotcombined',true,...
    'plotline',[[39;39;40;40;39],[-104;-103;-103;-104;-104]],...
    'plottabs',true,...
    'plotlimits',[0,0.05],...
    'plottitleadd',' (exampleAreaPlot)',...
    'returnonlydata',false,...
    'plotareacumulative',true,...
    'summarizearea',true);








