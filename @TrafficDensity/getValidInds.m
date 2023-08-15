function obj = getValidInds(obj)
% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
% Get valid indices for processing density data, based on input limitations
if isempty(obj.cell)
    warning('Must load data before processing; loading now');
    obj = LoadData(obj);
end
obj.validinds = true(size(obj.cell.Cell));
if numel(obj.timeofday) < 8 % Time of day - default is to process all 8 bins (assuming discretization is 3 hour bins)
    obj.validinds = obj.validinds.*ismember(obj.cell.Time,uint8(obj.timeofday));
end
if numel(obj.dayofweek) < 7 % Day of week - default is to process whole week
    wd = weekday(obj.cell.Date);
    obj.validinds = obj.validinds.*ismember(wd,obj.dayofweek);
end
if numel(obj.monthofyear) < 12 % Month - default is to process whole year
    [~,M,~,~,~,~] = datevec(double(obj.cell.Date));
    obj.validinds = obj.validinds.*ismember(M,obj.monthofyear);
end
if ~isempty(obj.ACcategory) % Aircraft category
    obj.validinds = obj.validinds.*ismember(obj.cell.Category,obj.ACcategory);
end
if ~isempty(obj.height) % Height/altitude
    obj.validinds = obj.validinds.*ismember(obj.cell.Height,obj.height);
end

% Reduce valid indices based on geographic extent of area or track
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
obj.cellLims = false(obj.jobdata.GRID_Y_NUM,obj.jobdata.GRID_X_NUM);
obj.cellLims(obj.cellLatLim(1):obj.cellLatLim(2),obj.cellLonLim(1):obj.cellLonLim(2)) = true;

% Get indices into original cell
indmat = reshape(1:obj.jobdata.GRID_Y_NUM*obj.jobdata.GRID_X_NUM,obj.jobdata.GRID_Y_NUM,obj.jobdata.GRID_X_NUM);
geoinds = indmat(obj.cellLatLim(1):obj.cellLatLim(2),obj.cellLonLim(1):obj.cellLonLim(2));
geoinds = geoinds(:);
georeduceinds = ismember(obj.cell.Cell,uint32(geoinds));
obj.validinds = obj.validinds.*georeduceinds;
end