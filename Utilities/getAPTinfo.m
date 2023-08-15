% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Purpose: get airport information and save for use by Traffic Density tool
% Airport information from FAA NASR: 
%    https://www.faa.gov/air_traffic/flight_info/aeronav/aero_data/NASR_Subscription/
% Process data from APT.txt record
% Data field location and description from Layout_Data/apt_rf.txt
%
function getAPTinfo(APTfilename)

% Download airport information if it does not exist
if ~exist('APTfilename','var')
    disp('Downloading airport information (APT.txt)...');
    datafolder = [getenv('TrafficDensityPath'),'/Data'];
    websave([datafolder, '/APT'],'https://nfdc.faa.gov/webContent/28DaySub/2020-10-08/APT.zip');
    unzip('APT.zip',datafolder);
    APTfilename = [datafolder, '/APT.txt'];
end

fid = fopen(APTfilename);

tic;
ii = 1;
clear ICAO lat lon el commercialAnnualOps
while true
    l = fgetl(fid);
    if l(1) == -1
        break;
    end
    
    % Check for airport    
    if ~strcmp(l(1:3),'APT') || ~strcmp(l(15:21),'AIRPORT') 
        continue;
    end

    ICAOstr = l(1211:1217);
    
    % Only keep airports with ICAO identifiers
    if all(ICAOstr==32)
        continue;
    end
    
    ICAO{ii} = ICAOstr(ICAOstr~=32); %#ok<AGROW> % Remove spaces and save
    
    commercialAnnualOps(ii) = str2double(l(1026:1031)); %#ok<AGROW> % Number of commerical ops could be used to reduce the number of airports
    
    % Location information
    latstr = l(539:550);
    lat(ii) = str2double(latstr(1:end-1))/3600; %#ok<AGROW>
    if latstr(end)=='S'
        lat(ii) = -lat(ii); %#ok<AGROW>
    end
    lonstr = l(566:577);
    lon(ii) = -str2double(lonstr(1:end-1))/3600; %#ok<AGROW>
    if lonstr(end)=='E' % Expect most to be 'W'
        lon(ii) = -lon(ii); %#ok<AGROW>
    end
    elstr = l(579:585);
    el(ii) = str2double(elstr); %#ok<AGROW>
    
    ii = ii+1;
end
toc;
fclose(fid);

% Save table output
APT = table(ICAO',lat',lon',el',commercialAnnualOps');
APT.Properties.VariableNames = {'ICAO','lat_deg','lon_deg','el_ft','commercialAnnualOps'};
save([getenv('TrafficDensityPath'),'/Data/APT.mat'],'APT');

