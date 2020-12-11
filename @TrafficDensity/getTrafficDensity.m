function obj = getTrafficDensity(obj)
% Copyright 2019 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Estimate the traffic density given indices (identified through user input)

% Define appropriate properties
jobdatastr = 'jobdata';
cellstr = 'cell';
termaplegendstr = 'termaplegend';

ncells = obj.(jobdatastr).GRID_X_NUM*obj.(jobdatastr).GRID_Y_NUM*obj.jobdata.GRID_H_NUM; % Total number of cells for 3D matrix
inds = obj.(cellstr).Cell(obj.validinds)+obj.(jobdatastr).GRID_X_NUM.*obj.(jobdatastr).GRID_Y_NUM.*uint32(obj.(cellstr).Height(obj.validinds)); % Indices for 3D matrix

% Get the number of observed hours
datetime = obj.(cellstr).Date(obj.validinds)*obj.jobdata.GRID_T_NUM+uint32(obj.(cellstr).Time(obj.validinds));
obj.obshrs = length(unique(datetime))/obj.jobdata.GRID_T_NUM*24;
clear datetime

% Compute density
% Get cell midpoints
obj.cellLatMidpoints = fliplr(obj.(termaplegendstr)(2)-1/obj.(termaplegendstr)(1)/2:-1/obj.(termaplegendstr)(1):obj.(termaplegendstr)(2)-obj.(jobdatastr).GRID_Y_NUM/obj.(termaplegendstr)(1));
obj.cellLonMidpoints = obj.(termaplegendstr)(3)+1/obj.(termaplegendstr)(1)/2:1/obj.(termaplegendstr)(1):obj.(termaplegendstr)(3)+obj.(jobdatastr).GRID_X_NUM/obj.(termaplegendstr)(1);

obj.cellLatCutpoints = fliplr(obj.(termaplegendstr)(2):-1/obj.(termaplegendstr)(1):obj.(termaplegendstr)(2)-obj.(jobdatastr).GRID_Y_NUM/obj.(termaplegendstr)(1));
obj.cellLonCutpoints = obj.(termaplegendstr)(3):1/obj.(termaplegendstr)(1):obj.(termaplegendstr)(3)+obj.(jobdatastr).GRID_X_NUM/obj.(termaplegendstr)(1);

% Get terrain
[lonm,latm] = meshgrid(obj.cellLonMidpoints,obj.cellLatMidpoints);
obj.el = ltln2val(obj.globe.Z,obj.globe.refvec,latm,lonm,'bilinear')*obj.f2m; % Convert to ft
obj.el(isnan(obj.el))=0; % Set to zero where no data exists

% Get cell area [NM^2]
[lon_offm,lat_offm] = meshgrid(obj.cellLonCutpoints(1:end-1),obj.cellLatCutpoints(1:end-1));
[lon_onm,lat_onm] = meshgrid(obj.cellLonCutpoints(2:end),obj.cellLatCutpoints(2:end));
cellarea = areaquad(lat_offm,lon_offm,lat_onm,lon_onm,obj.earthellipsoid); 

% Get cell height (including transition between AGL and MSL)
h_cutpoints = [obj.jobdata.AGL_LIMS,obj.jobdata.MSL_LIMS];
maxagl = find(h_cutpoints==max(obj.jobdata.AGL_LIMS));
dh = diff(h_cutpoints);
dh_mat = zeros(obj.(jobdatastr).GRID_Y_NUM, obj.(jobdatastr).GRID_X_NUM, obj.jobdata.GRID_H_NUM);
area_mat = zeros(obj.(jobdatastr).GRID_Y_NUM, obj.(jobdatastr).GRID_X_NUM, obj.jobdata.GRID_H_NUM);
for ii = 1:length(dh)
    area_mat(:,:,ii) = cellarea;
    if ii < maxagl % AGL
        dh_mat(:,:,ii) = dh(ii);
    elseif ii == maxagl % If the transition zone (from AGL to MSL) fix dh 
        dh_mat(:,:,ii) = dh(ii)-obj.el; % Any terrain elevation will reduce altitude range in transition zone       
    else % MSL only
        dh_mat(:,:,ii) = h_cutpoints(ii+1)-max(obj.el,h_cutpoints(ii)); % Limit minimum altitude cutpoint to be terrain elevation
    end
end
dh_mat(dh_mat<250) = nan; % Assume that need at least 250 ft altitude for reasonable density estimate

% Process each aircraft category separately (necessary when
% coupling with encounter model for collision rate estimate)
if obj.verbose
    % Only display message if processing 6 months or more of data
    if numel(obj.monthofyear)>=6
        disp('Computing density: computation may take some time, especially if computing for all data.'); 
    end
end
ACcategory = obj.ACcategory+1; % Make 1-based indexing for Matlab indexing

cat = obj.(cellstr).Category(obj.validinds);
TT = obj.(cellstr).TotalTime(obj.validinds);
MaxOcc = obj.(cellstr).MaxOcc(obj.validinds);
for ac = 1:length(ACcategory) 
    currACcat = ACcategory(ac); % 1 - discrete, 2 - VFR
    indsac = cat==currACcat-1;
    % Index into large arrays and tables requires nontrivial
    % time, so do it once here
    accuminds = inds(indsac);
    TotalTime = TT(indsac);
    % Aggregate results for each cell
    sumdensity = accumarray(accuminds,TotalTime,[ncells,1],@sum); % Observed aircraft seconds
    avgdensity = sumdensity/obj.obshrs/obj.hr2s; % Average number of aircraft in each cell    
    if obj.computemax % Run the mex file for computing the maximum density, which is considerably faster than Matlab built-in (see buildaccumarray.m under Utilities)
        maxdensity = double(accumarraymax_mex(accuminds,uint32(TotalTime),ncells))/24*obj.jobdata.GRID_T_NUM/obj.hr2s;
        maxocc = double(accumarraymax_mex(accuminds,uint32(MaxOcc(indsac)),ncells));
    
        % Reformat densities into matrices
        obj.count(currACcat).cmax = reshape(maxdensity,obj.(jobdatastr).GRID_Y_NUM,obj.(jobdatastr).GRID_X_NUM,obj.jobdata.GRID_H_NUM);
        obj.count(currACcat).cmaxocc = reshape(maxocc,obj.(jobdatastr).GRID_Y_NUM,obj.(jobdatastr).GRID_X_NUM,obj.jobdata.GRID_H_NUM);  
        
        % Following code is deprecated but included as example of how to
        % use Matlab built-in accumarray.m:
        % warning('Unable to run mex (compiled) function accumarraymax_mex; reverting to much slower Matlab built-in');
        % maxdensity = double(accumarray(accuminds,TotalTime,[ncells,1],@max))/24*obj.jobdata.GRID_T_NUM/obj.hr2s; % Maximum density per hour
        % maxocc = double(accumarray(accuminds,double(MaxOcc(indsac)),[ncells,1],@max)); % Maximum occupancy
    else
        if isfield(obj.count,'cmax')
            obj.count = rmfield(obj.count,'cmax');
        end
        if isfield(obj.count,'cmaxocc')
            obj.count = rmfield(obj.count,'cmaxocc');
        end
    end
    % Get upper bound of confidence interval - assumes poisson-distributed
    % data (see also Matlab built-in poissfit.m and private statpoisci.m)
    if obj.computeub
        lsum = obj.obshrs*obj.ciIndObsPerHr*avgdensity;
        k = lsum < 100;
        ub = zeros(size(avgdensity));
        if any(k); ub(k) = chi2inv(1-obj.cialpha/2, 2*(lsum(k)+1))/2; end
        k = ~k;    
        if any(k); ub(k) = norminv(1 - obj.cialpha/2,lsum(k),sqrt(lsum(k))); end
        ubavgdensity = ub./obj.obshrs/obj.ciIndObsPerHr;
        % Reformat density into matrix
        obj.count(currACcat).cub = reshape(ubavgdensity,obj.(jobdatastr).GRID_Y_NUM,obj.(jobdatastr).GRID_X_NUM,obj.jobdata.GRID_H_NUM);
    else
        if isfield(obj.count,'cub')
            obj.count = rmfield(obj.count,'cub');
        end
    end
    
    % Reformat densities into matrices
    obj.count(currACcat).c = reshape(avgdensity,obj.(jobdatastr).GRID_Y_NUM,obj.(jobdatastr).GRID_X_NUM,obj.jobdata.GRID_H_NUM);
    
    if obj.computestd % If want to evaluate the standard deviation, which requires more computation
        % Matlab is considerably faster at accumarray with @sum
        % than @std, so break computation of standard deviation
        % into separate accumarray @sum operations
        N = accumarray(accuminds,accuminds.*0+1,[ncells,1],@sum);
        den = N-1;
        den(den==0) = 1;
        den(den==-1) = 0;
        EX = sumdensity;
        EX2 = accumarray(accuminds,double(TotalTime).^2,[ncells,1],@sum);
        vardensity = (EX2-EX.^2./N)./den;
        vardensity(isnan(vardensity))=0;
        stddensity = (vardensity.^0.5)/24*obj.jobdata.GRID_T_NUM/obj.hr2s;
        %  stddensity = accumarray(inds(indsac),double(TotalTime)/24*obj.jobdata.GRID_T_NUM/obj.hr2s,[ncells,1],@std); % Numerically equivalent; left here as an example/benchmark
        obj.count(currACcat).cstd = reshape(stddensity,obj.(jobdatastr).GRID_Y_NUM,obj.(jobdatastr).GRID_X_NUM,obj.jobdata.GRID_H_NUM);        
    else
        if isfield(obj.count,'cstd')
            obj.count = rmfield(obj.count,'cstd');
        end
    end
    
    % Make values outside latitude/longitude limits NaN
    beyondLims = true([obj.(jobdatastr).GRID_Y_NUM,obj.(jobdatastr).GRID_X_NUM,obj.jobdata.GRID_H_NUM]);
    beyondLims(obj.cellLatLim(1):obj.cellLatLim(2),obj.cellLonLim(1):obj.cellLonLim(2),:) = false;
    for ff = fieldnames(obj.count(currACcat))'
        obj.count(currACcat).(ff{1})(beyondLims) = nan;
        if obj.correctcoverage % Correct for the cell coverage
            coveragecorrect = obj.cellcoverage;
            coveragecorrect(coveragecorrect<obj.noCoverageThreshold) = 0; % Set limit on correcting coverage
        else % Do not correct, but apply coverage threshold 
            coveragecorrect = double(obj.cellcoverage>=obj.noCoverageThreshold); % All 0's and 1's
        end
        coveragecorrect(coveragecorrect==0) = nan;
        obj.count(currACcat).(ff{1}) =  obj.count(currACcat).(ff{1})./coveragecorrect;
    end
    
    obj.density(currACcat).rho = obj.count(currACcat).c./(dh_mat.*obj.ft2nm.*area_mat); % AC Density (AC/NM^3)
    if obj.computemax 
        obj.density(currACcat).rho_max = obj.count(currACcat).cmax./(dh_mat.*obj.ft2nm.*area_mat); % Max average AC Density (AC/NM^3)
        obj.density(currACcat).rho_max_occ = obj.count(currACcat).cmaxocc./(dh_mat.*obj.ft2nm.*area_mat); % Max occupancy AC Density (AC/NM^3)
    else
        if isfield(obj.density,'rho_max')
            obj.density = rmfield(obj.density,'rho_max');
        end
        if isfield(obj.density,'rho_max_occ')
            obj.density = rmfield(obj.density,'rho_max_occ');
        end
    end
    if obj.computeub
        obj.density(currACcat).rhoub = obj.count(currACcat).cub./(dh_mat.*obj.ft2nm.*area_mat); % AC Density (AC/NM^3)
    else
        if isfield(obj.density,'rhoub')
            obj.density = rmfield(obj.density,'rhoub');
        end
    end
    if obj.computestd
        obj.density(currACcat).rho_std = obj.count(currACcat).cstd./(dh_mat.*obj.ft2nm.*area_mat); % AC Density (AC/NM^3)    
    else
        if isfield(obj.density,'rho_std')
            obj.density = rmfield(obj.density,'rho_std');
        end
    end
end

end