function loadSaveTerrain(datafolder,jobfile)
% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: to load and save terrain information (NOAA GLOBE)

% GLOBE tiles: https://www.ngdc.noaa.gov/mgg/topo/globeget.html
% ESRI headers: https://www.ngdc.noaa.gov/mgg/topo/elev/esri/hdr/

currDir = pwd;
cd(datafolder);
%% Download tiles and headers if they do not already exist
tiles = {'a10g','b10g','e10g','f10g'};
headers = {'a10g.hdr','b10g.hdr','e10g.hdr','f10g.hdr'};

% Download tiles and headers if they do not already exist
for i = 1:numel(tiles)
    if ~exist(tiles{i},'file')
        disp('Downloading GLOBE tiles...');
        if ispc
            websave(tiles{i},['https://www.ngdc.noaa.gov/mgg/topo/DATATILES/elev/' tiles{i} '.zip']);
            unzip([tiles{i} '.zip'],datafolder)
        else
            websave(tiles{i},['https://www.ngdc.noaa.gov/mgg/topo/DATATILES/elev/' tiles{i} '.gz']);
            untar([tiles{i} '.gz'],datafolder)
        end
    end

    if ~exist(headers{i},'file')
        disp('Downloading ESRI headers...');
        websave(headers{i},['https://www.ngdc.noaa.gov/mgg/topo/elev/esri/hdr/' headers{i}]);
    end
end

%% Process terrain information
if exist('jobfile','var') 
    td = TrafficDensity('jobfile',jobfile); 
else
    td = TrafficDensity;
end

% Get latitude/longitude limits of density data
latlim = [td.termaplegend(2)-td.jobdata.GRID_Y_NUM/td.termaplegend(1),td.termaplegend(2)];
lonlim = [td.termaplegend(3),td.termaplegend(3)+td.jobdata.GRID_X_NUM/td.termaplegend(1)];

% Matlab globedem function may error when latlim/lonlim are not integers
latlim(1) = floor(latlim(1));
latlim(2) = ceil(latlim(2));
lonlim(1) = floor(lonlim(1));
lonlim(2) = ceil(lonlim(2));

[Z,refvec] = globedem(datafolder,1,latlim,lonlim); % Read in elevation data

save([getenv('TrafficDensityPath'),'/Data/globe.mat'],'Z','refvec');
disp('Finished processing terrain data');
cd(currDir);