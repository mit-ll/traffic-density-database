function obj = runTrack(obj)
% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Process density and collision risk given a user specified track

% Provide guidance for update rate
difftime = diff(obj.track.Time_s(1:2));
if difftime<5 || difftime>30
    warning('Recommended time between updates is 5-30 s')
    if difftime<5
        warning('Too high an update rate will likely result in unnecessary computation')
    end
    if difftime>30
        warning('Too low an update rate may result in insufficeint fidelity')
    end
end

% Check that track has sufficient fields and equal sizes
tsize = size(obj.track.Time_s);
trackfields = fieldnames(obj.track);
for tt = 1:length(trackfields)
    if any(size(obj.track.(trackfields{tt}))~=tsize)
        error('All fields in track struction must have same size');
    end
end

h_cutpoints = [obj.jobdata.AGL_LIMS,obj.jobdata.MSL_LIMS];

jobdatastr = 'jobdata';

% Get latitude and longitude bins for each track point
[~,~,lati] = histcounts(obj.track.Latitude_deg,obj.cellLatCutpoints);
[~,~,loni] = histcounts(obj.track.Longitude_deg,obj.cellLonCutpoints);
tind2d = sub2ind([obj.(jobdatastr).GRID_Y_NUM,obj.(jobdatastr).GRID_X_NUM],lati,loni);

% Get correct altitude bins for each track point
elt = obj.el(tind2d);
altagl = obj.track.Altitude_MSL_ft-elt;
[~,~,alti] = histcounts(altagl,h_cutpoints); 
[~,~,altimsl] = histcounts(obj.track.Altitude_MSL_ft,h_cutpoints);
if ~isempty(obj.jobdata.MSL_LIMS) 
    useMSLalts = obj.track.Altitude_MSL_ft>=min(obj.jobdata.MSL_LIMS) & altagl>=max(obj.jobdata.AGL_LIMS); % If above AGL limits and in MSL region, use MSL bins
    alti(useMSLalts) = altimsl(useMSLalts);
    alti(altagl<min(obj.jobdata.AGL_LIMS)) = 1; % Below altitude limits
    alti(obj.track.Altitude_MSL_ft>max(obj.jobdata.MSL_LIMS)) = obj.jobdata.GRID_H_NUM; % Above altitude limits
else
    alti(altagl<min(obj.jobdata.AGL_LIMS)) = 1; % Below altitude limits
    alti(altagl>max(obj.jobdata.AGL_LIMS)) = obj.jobdata.GRID_H_NUM; % Above altitude limits
end

% Check for low-altitude trajectories that go above 5000 ft (threshold for
% change in the horizontal discretization)
if any(altagl<5000) && any(altagl>5000) && ~isempty(strfind(obj.filenames.cell,'lowh'))
    warning('Track goes above 5000 ft! lowh cell does not have coverage above 5000 ft. Traffic density estimates for track points above 5000 ft use data for the 3000-5000 ft altitude bin.');
end

% Get airspace class
obj = getTrackAirspace(obj); 

ACcategory = obj.ACcategory+1; % Make 1-based indexing for Matlab indexing
if isempty(ACcategory); ACcategory = [1,2]; end
obj.rate = [];
obj.relSpeed = [];
for ac = 1:length(ACcategory)
    currACcat = ACcategory(ac); % 1 - discrete, 2 - VFR/1200-code
    tind = sub2ind(size(obj.density(currACcat).rho),lati,loni,alti); % indices into density matix for each time
    if currACcat==1 % Get proper model characteristics property
        encModel = 'cor';
    elseif currACcat==2
        encModel = 'uncor';
    end
    
    % Get altitudes according to encounter model
    encModel_h_cutpoints = obj.(encModel).h_cutpoints;
    
    [~,~,alti_em] = histcounts(altagl,encModel_h_cutpoints);
    [~,~,altimsl_em] = histcounts(obj.track.Altitude_MSL_ft,encModel_h_cutpoints);
    useMSLalts_em = obj.track.Altitude_MSL_ft>=obj.(encModel).h_cutpoints(obj.(encModel).maxagl+1) & altagl>=obj.(encModel).h_cutpoints(obj.(encModel).maxagl);
    alti_em(useMSLalts_em) = altimsl_em(useMSLalts_em);
    alti_em(altagl<min(encModel_h_cutpoints)) = 1; % If lower than lower limit of encounter model, use lowest bin
	alti_em(obj.track.Altitude_MSL_ft>max(encModel_h_cutpoints)) = length(encModel_h_cutpoints)-1; 
    ualts_em = unique(alti_em);
    relSpeed = zeros(size(alti_em));
    for aa = 1:length(ualts_em) % For each unique altitude
        curralt = ualts_em(aa);
        currinds = find(alti_em==curralt);
        if curralt==0 % If out of bounds (too high for encounter model), continue
            relSpeed(currinds) = nan;
            continue;
        end
        speedcutpoints = obj.(encModel).boundaries{obj.(encModel).speedvar};
        p_speed_alt = obj.(encModel).p_speed_alt(:,curralt);
        
        nintrSpeed = 1000;
        intrSpeed = interp1([0;cumsum(p_speed_alt)],speedcutpoints,linspace(0.0001,0.9999,nintrSpeed)); % Get intruder speed deterministically (and avoid extreme outliers)
        ownSpeed = obj.track.Speed_kts(currinds);
        [uniqueSpeeds, ~, indexSpeeds] = unique(ownSpeed);
        ncurrinds = numel(uniqueSpeeds);
        if size(uniqueSpeeds,1)==1 % Make sure that is a column vector
            uniqueSpeeds = uniqueSpeeds';
        end
        intrSpeed_mat = repmat(intrSpeed,ncurrinds,1);
        ownSpeed_mat = repmat(uniqueSpeeds,1,nintrSpeed);
        
        relSpeedFun = @(heading)sqrt((intrSpeed_mat.*cos(heading)-ownSpeed_mat).^2+(intrSpeed_mat.*sin(heading)).^2);
        % SUGGESTION: consider speeding this up either with another
        % integration mechanism (e.g., sampling, trapz)
        relSpeedUnique = mean(integral(relSpeedFun,0,pi,'ArrayValued',true)/pi,2); % Average relative speed
        
        for bb = 1:ncurrinds
            relSpeed(currinds(indexSpeeds==bb))=relSpeedUnique(bb);
        end
    end
    obj.relSpeed{currACcat} = relSpeed;
        
    %Filter out NaN tind indices
    tind_notnan = tind(~isnan(tind)); %indices into density matix for each time that is not NaN
    ind_notnan = ~isnan(tind); %indices into each metric that is not NaN
    relSpeed = relSpeed(ind_notnan);
    
    % Get the collision rates
    obj.rate(currACcat).rateavg = NaN(size(tind));
    obj.rate(currACcat).rateavg(ind_notnan) = obj.density(currACcat).rho(tind_notnan).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;
    if obj.computemax
        obj.rate(currACcat).ratemax = NaN(size(tind));
        obj.rate(currACcat).ratemax(ind_notnan) = obj.density(currACcat).rho_max(tind_notnan).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;
        obj.rate(currACcat).ratemaxocc = NaN(size(tind));
        obj.rate(currACcat).ratemaxocc(ind_notnan) = obj.density(currACcat).rho_max_occ(tind_notnan).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;
    end
    if obj.computeub
        obj.rate(currACcat).rateavgub = NaN(size(tind));
        obj.rate(currACcat).rateavgub(ind_notnan) = obj.density(currACcat).rhoub(tind_notnan).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;
    end
    if obj.computestd
        obj.rate(currACcat).ratestd = NaN(size(tind));
        obj.rate(currACcat).ratestd(ind_notnan) = obj.density(currACcat).rho_std(tind_notnan).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;
    end
    
    % Update counts and density to correspond to track data
    obj.density(currACcat).rho = obj.density(currACcat).rho(tind_notnan);
    obj.count(currACcat).c = obj.count(currACcat).c(tind_notnan);
    if obj.computemax
        obj.density(currACcat).rho_max = obj.density(currACcat).rho_max(tind_notnan);
        obj.density(currACcat).rho_max_occ = obj.density(currACcat).rho_max_occ(tind_notnan);
        obj.count(currACcat).cmax = obj.count(currACcat).cmax(tind_notnan);
        obj.count(currACcat).cmaxocc = obj.count(currACcat).cmaxocc(tind_notnan);
    end
    if obj.computeub
        obj.density(currACcat).rhoub = obj.density(currACcat).rhoub(tind_notnan);
        obj.count(currACcat).cub = obj.count(currACcat).cub(tind_notnan);
    end
    if obj.computestd
        obj.density(currACcat).rho_std = obj.density(currACcat).rho_std(tind_notnan);
        obj.count(currACcat).cstd = obj.count(currACcat).cstd(tind_notnan);
    end
    
end

if obj.plotresults && length(intersect([1,2],ACcategory))==2
    obj.plot('plotcombined',true);
else
    obj.plot;
end

end

% Get airspace class for each track point
function obj = getTrackAirspace(obj)
    nairspace = numel(obj.airspace.CLASS);
   
    % Narrow down airspaces for testing based on intersection of bounding boxes
    airspaceBoundingBox = reshape([obj.airspace.BOUNDINGBOX_deg{:}],4,nairspace)';
    airspaceBoundingBoxPosVec = [airspaceBoundingBox(:,1),airspaceBoundingBox(:,2),airspaceBoundingBox(:,3)-airspaceBoundingBox(:,1),airspaceBoundingBox(:,4)-airspaceBoundingBox(:,2)];
    
    trackPosVec = [obj.lonlim(1),obj.latlim(1),diff(obj.lonlim),diff(obj.latlim)];
    
    airspaceIntersection = find(rectint(airspaceBoundingBoxPosVec,trackPosVec)>0);
    maxAirspaceAlt = max(obj.airspace.HIGHALT_ft_msl(airspaceIntersection));    
    
    possibleAirspace = table2cell(obj.airspace(airspaceIntersection,:));
    
    nairspaceInt = length(airspaceIntersection);
    ntrack = length(obj.track.Time_s);
    obj.airspaceClass = char(zeros(ntrack,1));
    
    % Get column numbers for cell
    varNames = obj.airspace.Properties.VariableNames;
    HIGHALT_ft_msl = strcmp(varNames,'HIGHALT_ft_msl');
    LOWALT_ft_msl = strcmp(varNames,'LOWALT_ft_msl');
    BOUNDINGBOX_deg = strcmp(varNames,'BOUNDINGBOX_deg');
    LON_deg = strcmp(varNames,'LON_deg');
    LAT_deg = strcmp(varNames,'LAT_deg');
    CLASS = strcmp(varNames,'CLASS');
    for tt = 1:ntrack
        currClass = [];
        currAlt = obj.track.Altitude_MSL_ft(tt);
        currLat = obj.track.Latitude_deg(tt);
        currLon = obj.track.Longitude_deg(tt);        
        if currAlt<=maxAirspaceAlt % Classify as B, C, or D if altitude less than maximum of all airspaces being checked
            for aa = 1:nairspaceInt % Check each possible airspace
                currAirspace = possibleAirspace(aa,:);
                % Check that altitude is within bounds
                if currAlt<=currAirspace{HIGHALT_ft_msl} && currAlt>=currAirspace{LOWALT_ft_msl} && ... % Altitudes
                        currLat<=currAirspace{BOUNDINGBOX_deg}(4) && currLat>=currAirspace{BOUNDINGBOX_deg}(2) && ... % Latitudes
                        currLon<=currAirspace{BOUNDINGBOX_deg}(3) && currLat>=currAirspace{BOUNDINGBOX_deg}(1) % Longitudes    
                        
                        % If desired, 3rd-party InPolygon can provide additonal speedups. See https://www.mathworks.com/matlabcentral/fileexchange/20754-fast-inpolygon-detection-mex.
                        inAirspace = inpolygon(currLon,currLat,currAirspace{LON_deg},currAirspace{LAT_deg}); % If in current airspace (use matlab build-in inpolygon.m)
                    if inAirspace 
                        currClass = currAirspace{CLASS};
                        break;
                    end
                end
            end
        end
        
        if isempty(currClass) % Classify as E, G, or A
            if currAlt >= 18000 
                currClass = 'A'; 
            else
                currClass = 'O'; % E or G
            end
        end
        obj.airspaceClass(tt) = currClass;
    end
end
