function [outdata,summarize] = plot(obj,varargin)
% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% PLOT method will display results of traffic density processing on a map
%   By default, the PLOT method will display the event rate information,
%   but can optionally plot the traffic density, observed traffic counts,
%   radar coverage, or airspace class. The latitude and longitude bounds
%   are set based on the user specified area or track TrafficDensity
%   properties.
%
%   Note: the following examples assume that an object `td` of type
%   TrafficDensity has been instantiated.
% 
%   td.PLOT(Name,Value) is the syntax for specifying the following options:
%    plottype - specify which data to plot from the following:
%       rate - event rate {default}
%       density - traffic density 
%       count - raw traffic counts from traffic database
%       coverage - estimated radar coverage
%       airspace - airspace class information
%           example: td.plot('plottype','density');
%    Note: see main README for a more detailed description.   
% 
%    plotscale - option for linear or logarithmic plotting of data:
%       log10 - base 10 logarithm {default}
%       linear - plot data as generated 
%           example: td.plot('plotscale','linear');
%
%    plotvalue - which summary statistic to display:
%       avg - average for each cell over time {default}
%       max - maximum for each cell over time
%       maxocc - maximum instantaneous occupancy
%       avgupperbound - the confidence interval upper bound for the average
%       std - standard deviation for the average (1 std above average is plotted)
%           example: td.plot('plotvalue','avgupperbound');
%    Note: see main README for a more detailed description.
%
%    plotcombined - whether to combine 1200-code/discrete and optionally
%    noncooperative results into a single additional aggregate estimate:
%       false - do not combine {default}
%       true - do combine with an additional set of figures
%            example: td.plot('plotcombined',true);
%     Note: plotcombined is simply a summation: thus, max and maxocc are 
%     not directly observed. 
%     Note: to include noncooperative plots, must specify TrafficDensity
%     class property `processNoncoop` as true
% 
%    plotline - option for overlaying a line on the map figures
%       If specified, must be an nx2 matrix with the columns corresponding
%       to latitude and longitude (in degrees):
%           example: td.plot('plotline',[[39;39;40;40;39],[-104;-103;-103;-104;-104]])
%
%    plottabs - option for specifying figures as tabs (to reduce the 
%    overall number of figures) or separate figures:
%       true - plot using tabs {default}
%       false - plot using only figures
%           example: td.plot('plottabs',false);
%
%    plotlimits - override the default calculation for determining minimum
%    and maximum limits for values; must be a 1x2 matrix
%           example: td.plot('plotlimits',[0,0.05]);
%
%    plottitleadd - option for specifying additional text to be appended to
%    the title string
%           example: td.plot('plottitleadd',' (Optional Additional Text)');
%
%    returnonlydata - will extract the necessary data as output without
%    plotting; this will select the data based on the latitude/longitude
%    and altitude limits specified by the user in the TrafficDensity
%    properties
%       true - extract and return data as output from method
%       false - continue to plot data in addition to returning data (if
%       output specified) {default}
%           example: [outdata,summarize] = td.plot('returnonlydata',true);
%       
%    plotareacumulative - if plotting a geographic area, will summarize the
%    data into cumulative densities for each airspace class and altitude
%    layer:
%       true - plot the cumulative density
%       false - do not plot the cumulative density {default}
%           example: td.plot('plotareacumulative',true);
%
%    summarizearea - whether to estimate aggregate summary statistics based
%    on altitude layer and airspace class if processing a geographic area:
%       true - summarize the area {default}
%       false - do not summarize the area
%    Note: will output data to screen if verbose TrafficDensity property
%    is true
%
%   [outdata,summarize] = td.PLOT(__) returns processed data:
%       outdata - extracted data as specified by user using the 'plottype'
%           and 'plotvalue' options; will conform to latitude, longitude,
%           and altitude limits specified by the user in the track and area
%           TrafficDensity properties
%       summarize - aggregated data based on airspace class and altitude
%
%   See also EXAMPLEAREAPLOT.

ncolors = 256; % Number of colors to plot

p = inputParser;
addParameter(p,'plottype','rate',@(x)ismember(x,{'rate','density','count','coverage','airspace'}));   % Which results/data to plot
addParameter(p,'plotscale','log10',@(x)ismember(x,{'linear','log10'}));         % Scale for the plotting
addParameter(p,'plotvalue','avg',@(x)ismember(x,{'avg','max','maxocc','avgupperbound','std'}));       % Which value to plot 
addParameter(p,'plotcombined',false,@(x)islogical(x));                          % Whether to plot combined 1200-code/discrete 
addParameter(p,'plotline',zeros(0,2),@(x)size(x,2)==2);                         % For specifying a line to plot over the map [lat,lon]
addParameter(p,'plottabs',true,@(x)islogical(x));                               % Whether to use tabs in figures (for fewer overall figures)
addParameter(p,'plotlimits',zeros(0,2),@(x)all(size(x)==[1,2]))                 % Plot limits for color range
addParameter(p,'plottitleadd','',@(x)ischar(x))                                 % User specified text to add to figure name
addParameter(p,'returnonlydata',false,@(x)islogical(x));                        % Whether to only return processed data and not plot (useful for computing noncooperative and aggregate)
addParameter(p,'plotareacumulative',false,@(x)islogical(x));                    % Whether to plot cumulative density
addParameter(p,'summarizearea',true,@(x)islogical(x));                          % Whether to summarize area and plot to command window
addParameter(p,'predominateThreshold',0,@(x)isnumeric(x) && isscalar(x) && (x >= 0) && (x<=1)); % Ignore cells where fraction of airspace is < threshold 
parse(p,varargin{:})

% Get altitude information
h_cutpoints = [obj.jobdata.AGL_LIMS,obj.jobdata.MSL_LIMS];
maxAGLind = length(obj.jobdata.AGL_LIMS);
h_types = repmat({'AGL'},length(h_cutpoints),1);
h_types(1:length(h_cutpoints)>maxAGLind) = {'MSL'};

evalheights = obj.height+1; % Convert to 1-based indexing used by Matlab
if isempty(evalheights)
    evalheights = 1:obj.jobdata.GRID_H_NUM;
end

% Get necessary data based on user selection
switch p.Results.plottype
    case 'rate' % Event rate
        data = obj.rate;
        if isempty(data); error('Must evaluate traffic density using `run` method before plotting results'); end
        avgField = 'rateavg';
        maxField = 'ratemax';
        maxoccField = 'ratemaxocc';
        avgUBField = 'rateavgub';
        stdField = 'ratestd';
        units = 'Events/Flight Hour';
        titleval = 'Event Rate';
    case 'density' % Traffic density
        data = obj.density;
        if isempty(data); error('Must evaluate traffic density using `run` method before plotting results'); end
        avgField = 'rho';
        maxField = 'rho_max';
        maxoccField = 'rho_max_occ';
        avgUBField = 'rhoub';
        stdField = 'rho_std';
        units = 'AC/NM^3';
        titleval = 'Aircraft Density';
    case 'count' % Traffic counts in each cell
        data = obj.count;
        if isempty(data); error('Must evaluate traffic density using `run` method before plotting results'); end
        avgField = 'c';
        maxField = 'cmax';
        maxoccField = 'cmaxocc';
        avgUBField = 'cub';
        stdField = 'cstd';
        units = 'AC';
        titleval = 'Aircraft Count';  
    case 'coverage' % Radar coverage information
        data = obj.cellcoverage;
        if isempty(data); error('Must load data into object before plotting coverage; either evaluate traffic density using `run` method or load data using `LoadData` method'); end
        cm = flipud(jet(ncolors));
        figtitle = ['Cell Radar Coverage',p.Results.plottitleadd];
        if p.Results.plottabs
            figure('Name',figtitle);
            tg = uitabgroup;
        end
        for hh=1:obj.jobdata.GRID_H_NUM
            tabtitle = sprintf('Altitudes: [%i ft %s,%i ft %s]', h_cutpoints(hh),h_types{hh},h_cutpoints(hh+1),h_types{hh+1});
            if p.Results.plottabs
                thistab = uitab(tg,'Title',tabtitle);
                axes('Parent',thistab);
            else
                figure('Name',figtitle);
                title(tabtitle);
            end
            cbarh = mapplot(obj,data(:,:,hh),cm,[0,1],p.Results.plotline,true);
            ylabel(cbarh,'Fraction of Cell with Radar Coverage');
        end
        return
    case 'airspace' % Airspace class information
        data = obj.cellAirspace;
        if isempty(data); error('Must load data into object before plotting airspace class information; either evaluate traffic density using `run` method or load data using `LoadData` method'); end
        cm = flipud(jet(ncolors));
        for aa=fieldnames(data)'
            figtitle = [sprintf('Airspace Class (%s)',aa{1}),p.Results.plottitleadd];
            if p.Results.plottabs
                figure('Name',figtitle);
                tg = uitabgroup;
            end
            for hh=1:obj.jobdata.GRID_H_NUM
                tabtitle = sprintf('Altitudes: [%i ft %s,%i ft %s]', h_cutpoints(hh),h_types{hh},h_cutpoints(hh+1),h_types{hh+1});
                if p.Results.plottabs
                    thistab = uitab(tg,'Title',tabtitle);
                    axes('Parent',thistab);
                else
                    figure('Name',figtitle);
                    sgtitle(tabtitle);
                end
                cbarh = mapplot(obj,data.(aa{1})(:,:,hh),cm,[0,1],p.Results.plotline,true);
                ylabel(cbarh,'Fraction of Cell with Airspace');
                title(aa{1})
            end
        end
        return
end

ACcategory = obj.ACcategory+1; % Make 1-based indexing for Matlab indexing
if isempty(ACcategory); ACcategory = [1,2]; end
if obj.processNoncoop; ACcategory = [ACcategory,3]; end
if p.Results.plotcombined && length(intersect([1,2],ACcategory))==2 % For plotcombined option
    ACcategory = [ACcategory,4]; 
elseif p.Results.plotcombined && length(intersect([1,2],ACcategory))~=2
    warning('Combined results cannot be plotted: must evaluate both types of aircraft category (ACcategory property) to combine results');
end 
nac = length(ACcategory); % Number of categories
throwCombineWarning = false;
datafieldnames = fieldnames(data);

summarize = struct;
outdata = [];
for ac=1:nac % For each type of aircraft    
    currACcat = ACcategory(ac); % 1 - discrete, 2 - VFR, 3 - noncoop, 4 - aggregate (discrete and 1200-code with or without noncoop)
    if currACcat==1 % Get proper model characteristics property
        categorytext = 'Discrete-Code Transponder';
    elseif currACcat==2
        categorytext = '1200-Code Transponder';
    elseif currACcat==3
        categorytext = 'Estimated Noncooperative';
        % Use 1200-code data, and apply noncoopFactor 
        for ff = 1:length(datafieldnames)
            data(currACcat).(datafieldnames{ff}) = data(2).(datafieldnames{ff})*obj.noncoopFactor;
        end
    elseif currACcat==4
        categorytext = 'Aggregate Aircraft (1200-code + Discrete)';        
        for ff = 1:length(datafieldnames)
            d1 = data(1).(datafieldnames{ff});
            d2 = data(2).(datafieldnames{ff});
            disnan = isnan(d1) & isnan(d2); % Preserve nan for plotting coverage vs. no observed data
            d1(isnan(d1)) = 0;
            d2(isnan(d2)) = 0;
            if obj.processNoncoop % Add noncooperatives if needed
                d2 = d2*(1+obj.noncoopFactor);
                categorytext = 'Aggregate Aircraft (1200-code + Discrete + Noncooperative)';    
            end
            data(currACcat).(datafieldnames{ff}) = d1+d2;
            data(currACcat).(datafieldnames{ff})(disnan) = nan;
        end
    end
    
    for ff = datafieldnames'
        if obj.processtrack % For track
            outdata(currACcat).(ff{1}) = data(currACcat).(ff{1}); %#ok<AGROW>
        else
            outdata(currACcat).(ff{1}) = data(currACcat).(ff{1})(obj.cellLatLim(1):obj.cellLatLim(2),obj.cellLonLim(1):obj.cellLonLim(2),evalheights); %#ok<AGROW>
        end
    end
    
    if p.Results.returnonlydata % If only want resulting data, then do not plot
        continue;
    end
    
    figtitle = [sprintf('%s for %s',titleval,categorytext),p.Results.plottitleadd];
    if p.Results.plottabs
        mapfh = figure('Name',figtitle);
        if ~obj.processtrack; tg = uitabgroup; end
        
        if p.Results.plotareacumulative && ~obj.processtrack
            cfh = figure('Name',[figtitle,' (Cumulative Density)']);
            tgc = uitabgroup; % Tabs for plotareacumulative
        end  
    end
    
    
    switch p.Results.plotvalue
        case 'avg'
            plotdata = data(currACcat).(avgField);
            plottitle = 'Average';            
        case 'max'
            if obj.computemax
                plotdata = data(currACcat).(maxField);
                plottitle = 'Maximum Average';
                throwCombineWarning = true;
            else
                error('Maximum was not computed: select different `plotvalue`');
            end
        case 'maxocc'
            if obj.computemax
                plotdata = data(currACcat).(maxoccField);
                plottitle = 'Maximum Occupancy';
                throwCombineWarning = true;
            else
                error('Maximum was not computed: select different `plotvalue`');
            end
        case 'avgupperbound'
            if obj.computeub
                plotdata = data(currACcat).(avgUBField);
                plottitle = 'Upper Bound for Average';
            else
                error('Upper bound was not computed: select different `plotvalue`');
            end
        case 'std'
            if obj.computestd
                plotdata = data(currACcat).(stdField);
                plottitle = 'Standard Deviation';
            else
                error('Standard deviation was not computed: select different `plotvalue`');
            end
    end
    
    if obj.processtrack % For track
        timeplot = (obj.track.Time_s-obj.track.Time_s(1))/obj.hr2s;
        % Plot map
        subplot('Position',[0.5703    0.1100    0.3347    0.8150])
        [~,APTh] = mapplot(obj);
        plotm(obj.track.Latitude_deg,obj.track.Longitude_deg,'k')        
        plotcurr = plotdata;
        noCovInds = isnan(plotcurr);
        noACInds = plotcurr==0;
        if strcmp(p.Results.plotscale,'log10')
             plotcurr = log10(plotcurr);              
        end
        plottimesdiff = timeplot(end)/10;
        plottimesdiff = ceil(plottimesdiff*12)/12; % round to 5 minute interval
        plottimes = 0:plottimesdiff:timeplot(end);
        plottimeslat = interp1(timeplot,obj.track.Latitude_deg,plottimes,'linear');
        plottimeslon = interp1(timeplot,obj.track.Longitude_deg,plottimes,'linear'); 
        sn = scatterm(obj.track.Latitude_deg(~noACInds & ~noCovInds),obj.track.Longitude_deg(~noACInds & ~noCovInds),30,plotcurr(~noACInds & ~noCovInds),'filled');  %#ok<NASGU>
        sz = scatterm(obj.track.Latitude_deg(noACInds),obj.track.Longitude_deg(noACInds),30,[0.5,0.5,0.5],'o','LineWidth',2);
        sc = scatterm(obj.track.Latitude_deg(noCovInds),obj.track.Longitude_deg(noCovInds),50,'k','x','LineWidth',2); 
        st = scatterm(plottimeslat,plottimeslon,100,[0.5,0,0.5],'d','LineWidth',2); 
        if ~isempty(plottimeslat) && ~isempty(plottimeslon) % Any airports in this region?
            ss = scatterm(plottimeslat(1),plottimeslon(1),200,[0.5,0,0.5],'s','LineWidth',2); 
        else
            ss = scatterm([],[]);
        end
        if ~isempty(p.Results.plotlimits) % Set color limits if needed
            minmaxplot = p.Results.plotlimits;
            if strcmp(p.Results.plotscale,'log10')
                minmaxplot = log10(minmaxplot);
            end
            set(gca,'CLim',minmaxplot);
        end
        legend([sz,sc,st,ss,APTh],'Zero with Coverage','No Coverage',sprintf('Every %i Minutes',plottimesdiff*60),'Track Start','Airports','Orientation','horizontal','Location','northoutside')
        legend('boxoff')
        cbarh = colorbar;
        title(plottitle)
        ylabel(cbarh,sprintf('%s (%s) [%s]',titleval,p.Results.plotscale,units))
        
        % Plot versus track time
        subplot('Position',[0.1300,0.3,0.3347,0.65]);    
        if obj.computeub
            if obj.computemax
                if obj.computestd
                    leg = semilogy(timeplot,data(currACcat).(avgField),'k',...
                        timeplot,data(currACcat).(avgUBField),'k-.',...
                        timeplot,data(currACcat).(maxField),'r',...
                        timeplot,data(currACcat).(maxoccField),'r-.',...
                        timeplot,data(currACcat).(avgField)+data(currACcat).(stdField),'b-.');
                else
                    leg = semilogy(timeplot,data(currACcat).(avgField),'k',...
                        timeplot,data(currACcat).(avgUBField),'k-.',...
                        timeplot,data(currACcat).(maxField),'r',...
                        timeplot,data(currACcat).(maxoccField),'r-.');
                end
            else
                if obj.computestd
                    leg = semilogy(timeplot,data(currACcat).(avgField),'k',...
                        timeplot,data(currACcat).(avgUBField),'k-.',...
                        timeplot,data(currACcat).(avgField)+data(currACcat).(stdField),'b-.');
                else
                    leg = semilogy(timeplot,data(currACcat).(avgField),'k',...
                        timeplot,data(currACcat).(avgUBField),'k-.');
                end
            end
        else
            if obj.computemax
                if obj.computestd
                    leg = semilogy(timeplot,data(currACcat).(avgField),'k',...
                        timeplot,data(currACcat).(maxField),'r',...
                        timeplot,data(currACcat).(maxoccField),'r-.',...
                        timeplot,data(currACcat).(avgField)+data(currACcat).(stdField),'b-.');
                else
                    leg = semilogy(timeplot,data(currACcat).(avgField),'k',...
                        timeplot,data(currACcat).(maxField),'r',...
                        timeplot,data(currACcat).(maxoccField),'r-.');
                end
            else
                if obj.computestd
                    leg = semilogy(timeplot,data(currACcat).(avgField),'k',...
                        timeplot,data(currACcat).(avgField)+data(currACcat).(stdField),'b-.');
                else
                    leg = semilogy(timeplot,data(currACcat).(avgField),'k');
                end
            end
        end
        xlabel('Track Time (hrs)')
        ylabel(sprintf('%s (%s) [%s]',titleval,p.Results.plotscale,units))
        xlim([min(timeplot),max(timeplot)])
        
        [leg, hasNoCoverage] = plotNoCoverage(plotcurr,timeplot,leg);
        set(gca,'Layer','top');
        labels = {'Average'};
        if obj.computeub
            labels{end+1} = 'Average Upper Bound'; %#ok<AGROW>
        end
        if obj.computemax
            labels{end+1} = 'Maximum 3 Hr Avg'; %#ok<AGROW> 
            labels{end+1} = 'Max Obs Occupancy'; %#ok<AGROW>
        end
        if obj.computestd
            labels{end+1} = '1 Std Dev above Avg'; %#ok<AGROW>
        end
        if hasNoCoverage
            labels{end+1} = 'No Radar Coverage'; %#ok<AGROW>
        end
        legend(leg,labels);
        text(0,min(ylim)*2,'Note: missing values indicate no observed data')
        grid on;
        
        subplot('Position',[0.1300,0.05,0.3347,0.15]);
        plotAirspaceClass(obj,timeplot,currACcat);
        if currACcat==1 || currACcat==4 % Discrete code intruder  
            xlabel('Note: for discrete-code intruder, gas-model only valid for VFR ownship in Class D, E, or G, and Discrete ownship in Class G.','FontSize',8);
        else % VFR 1200-code or noncooperative intruder
            xlabel('Note: for 1200-code/noncoop intruder, gas-model only valid for VFR ownship in Class C, D, E, or G, and Discrete ownship in Class D, E, or G.','FontSize',8);
        end        
        
    else % For area
        plotdataorig = plotdata;
        if isempty(p.Results.plotlimits)
            maxplot = max(plotdata(:)); % Get maximum over all altitudes for plotting
            minplot = min(plotdata(plotdata~=0));
        else
            maxplot = p.Results.plotlimits(2);
            minplot = p.Results.plotlimits(1);
            plotdata(plotdata>0 & plotdata<minplot & ~isnan(plotdata)) = minplot;
        end
        if strcmp(p.Results.plotscale,'log10')
            maxplot = log10(maxplot);
            minplot = log10(minplot);
        end
        
        plotnan = minplot-(maxplot-minplot)/(ncolors-2)*3/2; % No coverage
        plotnanlim = minplot-(maxplot-minplot)/(ncolors-2)*2; % Lower limit for plotting
        plotinf = minplot-(maxplot-minplot)/(ncolors-2)/2; % Zero density with coverage
        cm = jet(ncolors);
        cm(1,:) = [1,1,1];
        cm(2,:) = [0,0,0];       
        
        nalts = length(evalheights); % Altitudes to plot
        for hh=evalheights 
            % Plot map 
            minalt = h_cutpoints(hh);
            maxalt = h_cutpoints(hh+1);            
            
            tabtitle = sprintf('Altitudes: [%i ft %s,%i ft %s]',minalt,h_types{hh},maxalt,h_types{hh+1});
            if p.Results.plottabs
                figure(mapfh);
                thistab = uitab(tg,'Title',tabtitle);
                axes('Parent',thistab);
                tabtitle = [];
            else
                figure('Name',figtitle);
                title(tabtitle);
            end
            plotcurr = plotdata(:,:,hh);
            noCovInds = isnan(plotcurr);
            noACInds = plotcurr==0;
            if strcmp(p.Results.plotscale,'log10')
                plotcurr = log10(plotcurr);                
            end
            plotcurr(noCovInds) = plotnan;
            plotcurr(noACInds)= plotinf;            

            cbarh = mapplot(obj,plotcurr,cm,[plotnanlim,maxplot],p.Results.plotline);               
            ylabel(cbarh,sprintf('%s (%s) [%s]',titleval,p.Results.plotscale,units))
            if isempty(p.Results.plotlimits)
                whiteblack = 'White - No Coverage; Black - Zero with Coverage';
            else
                whiteblack = 'White - No Coverage; Black - Below Specified Limits with Coverage';
            end
            title({tabtitle,plottitle,whiteblack})
        end
        
        % Plot cumulative density if desired
        if p.Results.plotareacumulative   
            for aa = fieldnames(obj.cellAirspace)'
                cl = lines(nalts);
                
                if p.Results.plottabs
                    figure(cfh);
                    thistab = uitab(tgc,'Title',sprintf('Airspace Class: %s',aa{1}));
                    axes('Parent',thistab);
                else
                    figure('Name',figtitle);
                    title(tabtitle);
                end
                leg = {};
                for hh=evalheights
                    plotcurr = plotdataorig(:,:,hh);
                    valid = obj.cellAirspace.(aa{1})(:,:,hh) & ~isnan(plotcurr);
                    if sum(valid) < 3
                        continue;
                    end
                    w = obj.cellAirspace.(aa{1})(:,:,hh);
                    [histx,histy] = hist_w1d(plotcurr(valid),w(valid),'minval',0); % Get weighted distribution
                    leg = [leg;{sprintf('%i ft %s-%i ft %s',h_cutpoints(hh),h_types{hh},h_cutpoints(hh+1),h_types{hh+1})}]; %#ok<AGROW>
                    plot(histx,cumsum(histy),'Color',cl(hh,:)); hold on; % Plot cumulative distribution
                end
                grid on;
                ylim([0,1])
                xlabel(units)
                ylabel('Cumulative Density')
                legend(leg,'Location','southeast');
                
            end
        end
        
        % Summarize data if desired (output variable and to screen)
        if p.Results.summarizearea
            avgTable = table;
            for aa = fieldnames(obj.cellAirspace)'
                plotcurr = plotdataorig;
                as = obj.cellAirspace.(aa{1});
                
                valid = ~isnan(plotcurr);
                % Averages (weighted by airspace fraction for each cell)
                summarize(currACcat).avg.(aa{1}) =  squeeze(sum(plotcurr.*as.*valid,[1,2],'omitnan')./sum(as.*valid,[1,2],'omitnan'));
                
                % Maximums over regions with a given airspace                  
                % Estimate based on presence of an airspace in a cell
                summarize(currACcat).max.(aa{1}) =  squeeze(max(plotcurr.*(as.*valid>p.Results.predominateThreshold),[],[1,2],'omitnan'));

                avgTable.(aa{1}) = summarize(currACcat).avg.(aa{1});
            end
            
            % Get row names
            hrowname = cell(0);
            for hh = 1:length(h_cutpoints)-1
                hrowname = [hrowname;sprintf('%i ft %s,%i ft %s',h_cutpoints(hh),h_types{hh},h_cutpoints(hh+1),h_types{hh+1})]; %#ok<AGROW>
            end
            avgTable = avgTable(evalheights,:);
            avgTable.Properties.RowNames = hrowname(evalheights);
            
            if obj.verbose
                fprintf('Aggregate Average for %s %s (%s) [%s]:\n',titleval,plottitle,categorytext,units);  
                disp(avgTable);
            end
        end
    end
end

if throwCombineWarning && p.Results.plotcombined
    warning('Plotting both a maximum and ACcategory combined were specified; note that the combination is simply a summation, so the maximum combination is not observed in the data and may not be realistic');
end

end

function [cbarh,pAPTh] = mapplot(obj,data,cm,minmaxplot,plotline,plotdata)

plotDataOverride = false;
if exist('plotdata','var') && plotdata
    plotDataOverride = true;
end

if isempty(obj.cellLatLim) % If are running plot method before run method, determine the limits
    % Reduce valid indices based on geographic extent of area or
    % track
    if obj.processtrack % If processing a track
        obj.latlim = [min(obj.track.Latitude_deg),max(obj.track.Latitude_deg)];
        obj.lonlim = [min(obj.track.Longitude_deg),max(obj.track.Longitude_deg)];
    else % If processing an area
        obj.latlim = [min(obj.area.LatitudeLimit),max(obj.area.LatitudeLimit)];
        obj.lonlim = [min(obj.area.LongitudeLimit),max(obj.area.LongitudeLimit)];
    end
        
    % Get bin limits of area/track
    [~,~,obj.cellLatLim] = histcounts(obj.latlim,obj.cellLatCutpoints);
    [~,~,obj.cellLonLim] = histcounts(obj.lonlim,obj.cellLonCutpoints);
end

buf = 0.25; % Need a buffer due to Matlab not plotting all data otherwise
latlim = [obj.cellLatCutpoints(obj.cellLatLim(1)),obj.cellLatCutpoints(obj.cellLatLim(2)+1)];
lonlim = [obj.cellLonCutpoints(obj.cellLonLim(1)),obj.cellLonCutpoints(obj.cellLonLim(2)+1)];
latlim(1) = latlim(1)-buf; latlim(2) = latlim(2)+buf;
lonlim(1) = lonlim(1)-buf; lonlim(2) = lonlim(2)+buf;

termaplegend = [obj.termaplegend(1),obj.cellLatCutpoints(obj.cellLatLim(2)+1),obj.cellLonCutpoints(obj.cellLonLim(1))];

am = usamap(latlim,lonlim);
if ~exist('cm','var') || isempty(cm)
    cm = jet(256);
end
colormap(cm);
if ~obj.processtrack || plotDataOverride
    meshm(data(obj.cellLatLim(1):obj.cellLatLim(2),obj.cellLonLim(1):obj.cellLonLim(2)),termaplegend)
end
hold on
geoshow([obj.states.Lat],[obj.states.Lon],'Color',[0.8,0.8,0.8]);

if exist('minmaxplot','var') && ~isempty(minmaxplot) && ~any(isnan(minmaxplot)) && minmaxplot(1)~=minmaxplot(2)
    set(gca,'CLim',minmaxplot)
end

if ~obj.processtrack || plotDataOverride
    cbarh = colorbar;
else
    cbarh = [];
end

% Plot airports if desired
if obj.plotAirport
    validairport = obj.airport.lat_deg<latlim(2) & obj.airport.lat_deg>latlim(1) & obj.airport.lon_deg<lonlim(2) & obj.airport.lon_deg>lonlim(1);

    if sum(validairport)>obj.plotMaxNAirport
        validairport = find(validairport);
        [~,I] = sort(obj.airport.commercialAnnualOps(validairport),'descend','MissingPlacement','last');
        validairport = validairport(I(1:obj.plotMaxNAirport));
    end

    textm(obj.airport.lat_deg(validairport),obj.airport.lon_deg(validairport),obj.airport.ICAO(validairport),'FontSize',8)
    pAPTh = plotm(obj.airport.lat_deg(validairport),obj.airport.lon_deg(validairport),'k.');
else
    pAPTh = [];
end

if exist('plotline','var') && ~isempty(plotline) % Plot line if desired
    plotm(plotline(:,1),plotline(:,2),'Color',[0.8,0.8,0.8],'LineWidth',3)
end

set(get(am,'Title'), 'Visible', 'on')
tightmap('loose')
end


% Plot airspace class
function plotAirspaceClass(obj,timeplot,ac)
    AC = zeros(size(obj.airspaceClass));
    AC(obj.airspaceClass=='A') = 5;
    AC(obj.airspaceClass=='B') = 4;
    AC(obj.airspaceClass=='C') = 3;
    AC(obj.airspaceClass=='D') = 2;
    AC(obj.airspaceClass=='O') = 1;
    
    if ac==1 || ac==4 % Discrete code intruder
        r = AC==5 | AC==4 | AC==3;
        y = AC==2 | AC==1;
    else % 1200-code intruder
        r = AC==5 | AC==4;
        y = AC==3;
    end
    leg = {};
    
    rleg = plotPatch(r,timeplot,[],[0.5,5.5],[1,0.7,0.7]); 
    if ~isempty(rleg)
        leg = [leg,'Gas-Model Not Appropriate'];
    end
    yleg = plotPatch(y,timeplot,[],[0.5,5.5],[1,1,0.7]);
    if ~isempty(yleg)
        leg = [leg,'Gas-Model May Be Appropriate'];
    end
    
    p = plot(timeplot,AC,'k');
    ylim([0.5,5.5]);
    set(gca,'YTick',1:5);
    set(gca,'YTickLabel',{'E or G','D','C','B','A'});
    xlim([min(timeplot),max(timeplot)])
    set(gca,'XTickLabel',[])
    grid on;
    ylabel('Airspace Class')
    set(gca,'Layer','top','Box','on')
    leg = [leg,'Track Airspace Class'];
    legend([rleg,yleg,p],leg)
end

% Plot patches when signals are true
function leg = plotPatch(signal,xplot,leg,cylim,patchcolor)
    sigPad = [false;signal;false]; % Pad to identify no coverage at start or end
    startSig = find(diff(sigPad)>0);
    endSig = find(diff(sigPad)<0);
    nSig = length(startSig);
    for cc = 1:nSig
        s = max(startSig(cc)-1,1); % Correct start time due to padding
        e = min(endSig(cc),length(xplot));
        x = xplot([s,e,e,s]);
        y = [cylim(1),cylim(1),cylim(2),cylim(2)];
        p = patch(x,y,patchcolor,'EdgeColor','none'); hold on;
        if cc==1
            leg = [leg;p]; %#ok<AGROW>
        end
    end
end

% Determine and plot times with no coverage
function [leg, hasNoCoverage] = plotNoCoverage(coveragecorrect,timeplot,leg)    
    noCovInds = isnan(coveragecorrect); % No coverage
    cylim = ylim;
    leg = plotPatch(noCovInds,timeplot,leg,cylim,[0.8,0.8,0.8]);
    hasNoCoverage = sum(noCovInds)~=0;
end
