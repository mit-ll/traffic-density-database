function obj = SummarizeData(obj)
% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Summarize the loaded data, including data span and aggregate map
toplot = length(dbstack)==1; % Determine if called from command line (rather than another function/method) for plotting
if isempty(obj.cell)
    warning('Must load data before summarizing; loading now');
    obj = LoadData(obj);
end

ncells = obj.jobdata.GRID_Y_NUM*obj.jobdata.GRID_X_NUM;

udates = unique(obj.cell.Date); % Get unique dates
ndates = max(udates)-min(udates)+1; % Get theoretical date span

% Get cells per day
flightHoursPerDay = accumarray(obj.cell.Date-min(udates)+1,obj.cell.TotalTime,[ndates,1]);

% Plot
if toplot
    figure;
    plot(1:ndates,flightHoursPerDay/obj.hr2s,'.')
    title('Flight Hours Per Day')
    xlabel('Day from Start')
    ylabel('Flight Hours Per Day');
    grid on;
end

% Plot flight hours on map
flightHoursPerCell = accumarray(obj.cell.Cell,obj.cell.TotalTime,[ncells,1]);
if toplot
    m = zeros(obj.jobdata.GRID_Y_NUM, obj.jobdata.GRID_X_NUM);
    
    m(1:ncells) = flightHoursPerCell ./ obj.hr2s; % convert to hours
    m(m == 0) = nan;
    
    m(~isnan(m)) = log10(m(~isnan(m)));
    
    states=shaperead('usastatelo','UseGeoCoords',true, ...
        'Selector', ...
        {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})),'Name'});
    
    figure;
    usamap(m, obj.termaplegend);
    clmo surface
    cm = jet(256);
    cm(1,:) = [1,1,1];
    colormap(cm);
    meshm(m,obj.termaplegend)
    hold on
    h = colorbar;
    geoshow([states.Lat],[states.Lon],'Color','white');
    set(gca,'CLim',[-1,4])
    htick = get(h,'YTick');
    set(h,'YTickLabel',num2cell(10.^htick'),'YTick',htick);
    ylabel(h,'Flight Hours')
    title('Total Flight Time Observed')
end
end