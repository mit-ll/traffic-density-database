% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% getCellAirspace will compute airspaces for each cell given TrafficDensity
% object

td = TrafficDensity; % Instantiate object to get jobdata and airspace filename

% Load necessary airspace data files
tmp = load(td.filenames.airspace);
airspace = tmp.airspace;
airspaceCell = table2cell(airspace); % Use cell for quicker indexing
nairspace = numel(airspace.CLASS);

if any(airspace.CLASS=='E')
    error('getCellAirspace does not support Class E airspace');
end

% Get grid cell information
ncells2d = td.jobdata.GRID_X_NUM*td.jobdata.GRID_Y_NUM; % Number of 2d cells
ncells = td.jobdata.GRID_X_NUM*td.jobdata.GRID_Y_NUM*td.jobdata.GRID_H_NUM; % Total number of cells for 3D matrix

cellLatMidpoints = fliplr(td.termaplegend(2)-1/td.termaplegend(1)/2:-1/td.termaplegend(1):td.termaplegend(2)-td.jobdata.GRID_Y_NUM/td.termaplegend(1));
cellLonMidpoints = td.termaplegend(3)+1/td.termaplegend(1)/2:1/td.termaplegend(1):td.termaplegend(3)+td.jobdata.GRID_X_NUM/td.termaplegend(1);

cellLatCutpoints = fliplr(td.termaplegend(2):-1/td.termaplegend(1):td.termaplegend(2)-td.jobdata.GRID_Y_NUM/td.termaplegend(1));
cellLonCutpoints = td.termaplegend(3):1/td.termaplegend(1):td.termaplegend(3)+td.jobdata.GRID_X_NUM/td.termaplegend(1);

h_cutpoints = [td.jobdata.AGL_LIMS,td.jobdata.MSL_LIMS];
maxagl = find(h_cutpoints==max(td.jobdata.AGL_LIMS));

% Get terrain
ft2m = 3.28084; % Conversion factor (1 m = 3.28084 ft)
tmp = load('globe.mat');
globe.Z = tmp.Z;
globe.refvec = tmp.refvec;

[lonm,latm] = meshgrid(cellLonMidpoints,cellLatMidpoints);
el = ltln2val(globe.Z,globe.refvec,latm,lonm,'bilinear')*ft2m; % Convert to ft
el(isnan(el))=0; % Set to zero where no data exists

% Get column numbers when convert to cell
varNames = airspace.Properties.VariableNames;
HIGHALT_ft_msl = strcmp(varNames,'HIGHALT_ft_msl');
LOWALT_ft_msl = strcmp(varNames,'LOWALT_ft_msl');
HIGHALT_ft_agl = strcmp(varNames,'HIGHALT_ft_agl');
LOWALT_ft_agl = strcmp(varNames,'LOWALT_ft_agl');
BOUNDINGBOX_deg = strcmp(varNames,'BOUNDINGBOX_deg');
LON_deg = strcmp(varNames,'LON_deg');
LAT_deg = strcmp(varNames,'LAT_deg');
CLASS = strcmp(varNames,'CLASS');

% Get airspace bounding boxes (for initially reducing number of airspaces
% to check overlap with cells)
airspaceBoundingBox = reshape([airspace.BOUNDINGBOX_deg{:}],4,nairspace)';
airspaceBoundingBoxPosVec = [airspaceBoundingBox(:,1),airspaceBoundingBox(:,2),airspaceBoundingBox(:,3)-airspaceBoundingBox(:,1),airspaceBoundingBox(:,4)-airspaceBoundingBox(:,2)];

%% For each grid cell, get airspace information
tic;
clear airspaceClasses;
airspaceClasses = {'A','B','C','D','O'};
for aa = 1:length(airspaceClasses)
    cellAirspace.(airspaceClasses{aa}) = zeros([td.jobdata.GRID_Y_NUM,td.jobdata.GRID_X_NUM,td.jobdata.GRID_H_NUM]);
end
for cc = 1:ncells2d
    
    % Get cell bounding box
    [latsub,lonsub] = ind2sub([td.jobdata.GRID_Y_NUM,td.jobdata.GRID_X_NUM],cc);
    cellBoundingBox = [cellLonCutpoints(lonsub),cellLatCutpoints(latsub),cellLonCutpoints(lonsub+1),cellLatCutpoints(latsub+1)];
    cellBoundingBoxPosVec = [cellBoundingBox(:,1),cellBoundingBox(:,2),cellBoundingBox(:,3)-cellBoundingBox(:,1),cellBoundingBox(:,4)-cellBoundingBox(:,2)];
    cellel = el(latsub,lonsub);
    
    % Downsample number of airspaces based on cell
    airspaceIntersection = find(rectint(airspaceBoundingBoxPosVec,cellBoundingBoxPosVec)>0);    
    nairspaceInt = length(airspaceIntersection);
    possibleAirspace = airspaceCell(airspaceIntersection,:);
    
    % Loop over each altitude cell
    for hh = 1:td.jobdata.GRID_H_NUM
        hCellmin = h_cutpoints(hh);
        hCellmax = h_cutpoints(hh+1);
        
        if hCellmin<18000 % Assume that there is an altitude cutpoint at 18000 ft MSL (Class A boundary)
            cellAirspace.O(latsub,lonsub,hh) = 1;
        else
            cellAirspace.A(latsub,lonsub,hh) = 1;
        end
        
        % If cell defined by AGL, convert to MSL
        if hh<=maxagl 
            hCellmin = hCellmin+cellel;
        end
        if hh+1<=maxagl 
            hCellmax = hCellmax+cellel;
        end
        
        if hCellmax<=hCellmin || nairspaceInt==0 % If bin has no height or no overlapping airspaces, go to next bin (do nairspaceInt check here to set O and A airspaces first)
            continue;
        end
                
        for aa = 1:nairspaceInt % Check each possible airspace

            currAirspace = possibleAirspace(aa,:);         
                        
            % Check altitude
            hAirspacemin = currAirspace{LOWALT_ft_msl};
            hAirspacemax = currAirspace{HIGHALT_ft_msl};
            
            hOverlapmax = min(hAirspacemax,hCellmax);
            hOverlapmin = max(hAirspacemin,hCellmin);
            
            if hOverlapmax<=hOverlapmin % If there is no altitude overlap
                continue;
            end
            hfraction = (hOverlapmax-hOverlapmin)/(hCellmax-hCellmin); % Altitude fraction (to be combined with horizontal fraction)
            
            cellpoly = polyshape([cellBoundingBox(1),cellBoundingBox(3),cellBoundingBox(3),cellBoundingBox(1)],...
                [cellBoundingBox(2),cellBoundingBox(2),cellBoundingBox(4),cellBoundingBox(4)]);
            airspacepoly = polyshape(currAirspace{LON_deg},currAirspace{LAT_deg});
            polyout = intersect(cellpoly,airspacepoly);
            polyintersect = area(polyout);
            hozfraction = polyintersect/(cellBoundingBoxPosVec(3)*cellBoundingBoxPosVec(4));
            
            fraction = hozfraction.*hfraction;
            currClass = char(currAirspace{CLASS});
            
            cellAirspace.(currClass)(latsub,lonsub,hh)= cellAirspace.(currClass)(latsub,lonsub,hh)+fraction; % Assumes that airspaces do not overlap (mutually exclusive)
            cellAirspace.(currClass)(latsub,lonsub,hh) = min(cellAirspace.(currClass)(latsub,lonsub,hh),1); % Limit airspace fraction in cases where greater than 1 (in case of overlapping airspaces)
        end
        
        % Address Class O (E and G) and A
        if hCellmin<18000 % Assume that there is an altitude cutpoint at 18000 ft
            cellAirspace.O(latsub,lonsub,hh) = max(0,1-cellAirspace.B(latsub,lonsub,hh)-cellAirspace.C(latsub,lonsub,hh)-cellAirspace.D(latsub,lonsub,hh));
        else
            cellAirspace.A(latsub,lonsub,hh) = 1;
        end
    end
    fprintf('%i/%i\n',cc,ncells2d)
end

% Ensure that sum over classes for each spatial bin is 1
% Can be greater than 1 in cases where there is overlapping different 
% classes (B, C, D)
% Get sum over airspace classes
airspaceSum = zeros(size(cellAirspace.A));
for aa = 1:length(airspaceClasses)
    airspaceSum = airspaceSum+cellAirspace.(airspaceClasses{aa});
end
% Correct by normalizing
for aa = 1:length(airspaceClasses)
    cellAirspace.(airspaceClasses{aa}) = cellAirspace.(airspaceClasses{aa})./airspaceSum;
end
toc;

% Save results
save([getenv('TrafficDensityPath'),'/Data/cellAirspace.mat'],'cellAirspace');

% Plot results for verification 
figure('Name','Airspace Class Cell Verification (Cell Fraction for each Airspace)'); 
tg = uitabgroup;
ncol = floor(sqrt(td.jobdata.GRID_H_NUM));
nrow = ceil(td.jobdata.GRID_H_NUM/ncol);
cm = jet(256);
cm(1,:) = [1,1,1];
colormap(cm);
for aa = 1:length(airspaceClasses)
    thistab = uitab(tg,'Title',sprintf('Airspace Class: %s',airspaceClasses{aa}));
    axes('Parent',thistab);
    for hh = 1:td.jobdata.GRID_H_NUM
        hCellmin = h_cutpoints(hh);
        hCellmax = h_cutpoints(hh+1);
        subplot(nrow,ncol,hh)
        imagesc(cellLonMidpoints,cellLatMidpoints,squeeze(cellAirspace.(airspaceClasses{aa})(:,:,hh))); 
        
        set(gca,'YDir','normal');
        xlabel('Longtitude (deg)')
        ylabel('Latitude (deg)')
        title(sprintf('Altitudes: [%i,%i]',hCellmin,hCellmax));
        set(gca,'CLim',[0,1])
        axis image;
        grid on;
        colorbar;
    end
end
