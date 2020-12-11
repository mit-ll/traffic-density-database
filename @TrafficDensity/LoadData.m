function obj = LoadData(obj)
% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Load density, coverage, airspace, and terrain data
%
% Do not automatically load data when constructing object so that
% properties can be set and explored without the time to load the data

% Return if data has been previously loaded
if ~isempty(obj.cell)
    return;
end

%% Check whether user has sufficient memory 
[userview, systemView] = memory; % Get available memory
reqMemory_GB = 5.5; % 5.26 GB is required load data table (round to 5.5 GB)
reqRAM_GB = 18; % Up to 18GB RAM is required for processing
if userview.MaxPossibleArrayBytes/1e9 <= reqMemory_GB
    warning('At least %i GB of memory is required to load the traffic density tables', reqMemory_GB);
end
if systemView.PhysicalMemory.Total/1e9 <= reqRAM_GB
    warning('At least %i GB of memory is required to run the traffic density tool', reqRAM_GB);
end

a = tic;
if obj.verbose; disp('Loading traffic density data (this may take some time)'); end

%% Load density data
if ~exist(obj.filenames.cell,'file')
    error('Cell (traffic density) file does not exist on Matlab path');
end
tmp = load(obj.filenames.cell,'cell');
if isempty(strfind(obj.filenames.cell,'lowh'))
    if any(size(tmp.cell)~=[337086522, 7]) % Check data has expected size
        warning('Size of traffic density cell is different than expected (Expected: 301017469x7 table)');
    end
else
    if any(size(tmp.cell)~=[410100649, 7]) % Check data has expected size (low-altitude table)
        warning('Size of traffic density cell is different than expected (Expected: 410100649x7 table)');
    end
end
obj.cell = tmp.cell;

%% Load radar coverage data
if ~exist(obj.filenames.cellCoverage,'file')
    error('Radar coverage file does not exist on Matlab path');
end
tmp = load(obj.filenames.cellCoverage,'cellcoverage');
if isempty(strfind(obj.filenames.cell,'lowh'))
    if any(size(tmp.cellcoverage)~=[163, 373, 9]) % Check data has expected size
        warning('Size of radar coverage data is different than expected (Expected: 163x373x9)');
    end
else
    if any(size(tmp.cellcoverage)~=[652, 1492, 4]) % Check data has expected size (low-altitude table)
        warning('Size of radar coverage data is different than expected (Expected: 652x1492x4)');
    end
end
obj.cellcoverage = tmp.cellcoverage;

%% Load NOAA GLOBE terrain data
tmp = load('globe.mat');
if any(size(tmp.Z)~=[3481, 7681]) || numel(tmp.refvec)~=3 % Check data has expected size
    warning('Size of GLOBE terrain is different than expected (Expected Z: 3481x7681, Expected refvec: 1x3)');
end
obj.globe = tmp;

%% Load encounter model characteristics
obj = obj.loadEncModel;

%% Load airport information
tmp = load(obj.filenames.airport);
if any(size(tmp.APT)~=[2643, 5]) % Check data has expected size
    warning('Size of airport data is different than expected (Expected: 2643x5)');
end
obj.airport = tmp.APT;

%% Load airspace information
tmp = load(obj.filenames.airspace);
if any(size(tmp.airspace)~=[1249, 11]) % Check data has expected size
    warning('Size of airspace information is different than expected (Expected: 1249x11)');
end
obj.airspace = tmp.airspace;
tmp = load(obj.filenames.cellAirspace);
if isempty(strfind(obj.filenames.cellAirspace,'lowh'))
    if any(size(tmp.cellAirspace.A)~=[163, 373, 9]) || any(size(tmp.cellAirspace.B)~=[163, 373, 9]) ...
            || any(size(tmp.cellAirspace.C)~=[163, 373, 9]) || any(size(tmp.cellAirspace.D)~=[163, 373, 9])...
            || any(size(tmp.cellAirspace.O)~=[163, 373, 9]) % Check data has expected size
        warning('Size of cell airspace information is different than expected (Expected for all airspaces: 163x373x9)');
    end
else
    if any(size(tmp.cellAirspace.A)~=[652, 1492, 4]) || any(size(tmp.cellAirspace.B)~=[652, 1492, 4]) ...
            || any(size(tmp.cellAirspace.C)~=[652, 1492, 4]) || any(size(tmp.cellAirspace.D)~=[652, 1492, 4])...
            || any(size(tmp.cellAirspace.O)~=[652, 1492, 4]) % Check data has expected size (low-altitude table)
        warning('Size of cell airspace information is different than expected (Expected for all airspaces: 652x1492x4)');
    end
end
obj.cellAirspace = tmp.cellAirspace;

if obj.verbose
    GetSize(obj);
    fprintf('Time required to load data: %.2f s\n',toc(a));
end

%% Load mapping states data
obj.states = shaperead('usastatelo','UseGeoCoords',true,'Selector',{@(name) ~any(strcmp(name,{'Alaska','Hawaii'})),'Name'});

end
