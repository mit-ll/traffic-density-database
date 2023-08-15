function obj = runArea(obj)
% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Process density and collision risk given a user specified area

% Get heights to process
evalheights = obj.height+1; % Convert to 1-based indexing used by Matlab
if isempty(evalheights)
    evalheights = 1:obj.jobdata.GRID_H_NUM;
end
h_cutpoints = [obj.jobdata.AGL_LIMS,obj.jobdata.MSL_LIMS];

ACcategory = obj.ACcategory+1; % Make 1-based indexing for Matlab indexing
if isempty(ACcategory); ACcategory = [1,2]; end
obj.rate = [];
obj.summarize = [];
for ac = 1:length(ACcategory) % For each AC category being considered
    currACcat = ACcategory(ac); % 1 - discrete, 2 - VFR
    if currACcat==1 % Get proper model characteristics property
        encModel = 'cor';
    elseif currACcat==2
        encModel = 'uncor';
    end
    encModel_h_cutpoints = obj.(encModel).h_cutpoints;
    obj.rate(currACcat).rateavg = nan(size(obj.density(currACcat).rho));
    obj.rate(currACcat).ratemax = nan(size(obj.density(currACcat).rho_max));
    obj.rate(currACcat).ratemaxocc = nan(size(obj.density(currACcat).rho_max_occ));
    if obj.computeub
        obj.rate(currACcat).rateavgub = nan(size(obj.density(currACcat).rhoub));
    end
    if obj.computestd
        obj.rate(currACcat).ratestd = nan(size(obj.density(currACcat).rho_std));
    end
    for hh = evalheights % For each height in the density data to evaluate
        maxalt = h_cutpoints(hh+1);
        % Compute the relative speed
        % Get encounter model bin corresponding to traffic
        % density altitude bin: assume worst case (highest)
        % SUGGESTION: consider weighting the encounter model altitudes
        % rather than selecting the worst case
        worstCaseEncModelBin = find(encModel_h_cutpoints<maxalt,1,'last');
        if isempty(worstCaseEncModelBin)
            worstCaseEncModelBin = 1;
        end
        if worstCaseEncModelBin>length(encModel_h_cutpoints)-1 % If out of bounds of the encounter model, continue
            continue;
        end
        
        speedcutpoints = obj.(encModel).boundaries{obj.(encModel).speedvar};
        p_speed_alt = obj.(encModel).p_speed_alt(:,worstCaseEncModelBin);
        
        intrSpeed = interp1([0;cumsum(p_speed_alt)],speedcutpoints,linspace(0.0001,0.9999,10000)); % Get intruder speed deterministically (and avoid extreme outliers)
        if isempty(obj.ownspeed) % If ownspeed undefined, set at mean
            ownSpeed = mean(intrSpeed);
        elseif numel(obj.ownspeed)==1 % If own speed defined as scalar
            ownSpeed = obj.ownspeed;
        elseif numel(obj.ownspeed)==obj.jobdata.GRID_H_NUM % If own speed defined by altitudes
            ownSpeed = obj.ownspeed(hh);
        else
            error('For area processing, the ownspeed property must be defined as a scalar or must have the same number of elements as the number of altitude bins');
        end
        relSpeedFun = @(heading)sqrt((intrSpeed.*cos(heading)-ownSpeed).^2+(intrSpeed.*sin(heading)).^2);
        relSpeed = mean(integral(relSpeedFun,0,pi,'ArrayValued',true)/pi); % Average relative speed
        
        % Get the collision rates   
        obj.rate(currACcat).rateavg(:,:,hh) = obj.density(currACcat).rho(:,:,hh).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;
        obj.rate(currACcat).ratemax(:,:,hh) = obj.density(currACcat).rho_max(:,:,hh).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;
        obj.rate(currACcat).ratemaxocc(:,:,hh) = obj.density(currACcat).rho_max_occ(:,:,hh).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2; 
        if obj.computeub
            obj.rate(currACcat).rateavgub(:,:,hh) = obj.density(currACcat).rhoub(:,:,hh).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;     
        end
        if obj.computestd
            obj.rate(currACcat).ratestd(:,:,hh) = obj.density(currACcat).rho_std(:,:,hh).*4.*obj.macR.*obj.macH.*relSpeed*obj.ft2nm^2;
        end
    end        
end

if obj.plotresults
    obj.plot('plotcombined',true,'summarizearea',obj.verbose);
end

