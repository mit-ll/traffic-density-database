function generateTestData
% Copyright 2019 - 2023, MIT Lincoln Laboratory
% SPDX-License-Identifier: X11
%
%The purpose of this function is to generate unit test data for the
%TrafficDensity tool by manipulating the values of TotalTime and MaxOcc in
%the data created in generateTestCellData.m. In the test data, TotalTime
%and MaxOcc have the same values.

% Load in test data generated in generateTestCellData.m
tmp = load('cellOrigTest.mat','cell');
cell = tmp.cell;
origCell = cell;

saveDirectory = [getenv('TrafficDensityPath') filesep 'Testing' filesep 'TestData' filesep]; %Directory where plots will be saved

%% Default Test Data
%Create simple default test dataset where the density is solely dependent
%on altitude layer
testData = origCell;

% Redefine the TotalTime category based on height:
% Bin 0 corresponds to 1
% Bin 1 corresponds to 2
% Bin 2 corresponds to 3
testData.TotalTime(testData.Height==0) = uint32(1);
testData.TotalTime(testData.Height==1) = uint32(2);
testData.TotalTime(testData.Height==2) = uint32(3);
testData.MaxOcc = testData.TotalTime;
cell = testData; 

save([saveDirectory, 'testData.mat'],'cell','-v7.3');
clear cell;

%% Airspace Space Class Test Data
% Test data is 10x10 with 3 altitude bins
airspace.A = zeros(10,10,3);
airspace.B = zeros(10,10,3);
airspace.C = zeros(10,10,3);
airspace.D = zeros(10,10,3);
airspace.O = zeros(10,10,3);

%Data where all airspaces are equally represented in each cell
airspaceAll = airspace;
airspaceAll.A = ones(size(airspaceAll.A))*0.2;
airspaceAll.B = ones(size(airspaceAll.B))*0.2;
airspaceAll.C = ones(size(airspaceAll.C))*0.2;
airspaceAll.D = ones(size(airspaceAll.D))*0.2;
airspaceAll.O = ones(size(airspaceAll.O))*0.2;
cellAirspace = airspaceAll;
save([saveDirectory,'testAirspaceAll.mat'],'cellAirspace','-v7.3');

%Data with only class A airspace
airspaceA = airspace;
airspaceA.A = ones(size(airspaceA.A));
airspaceA.B = zeros(size(airspaceA.B));
airspaceA.C = zeros(size(airspaceA.C));
airspaceA.D = zeros(size(airspaceA.D));
airspaceA.O = zeros(size(airspaceA.O));
cellAirspace = airspaceA;
save([saveDirectory,'testAirspaceA.mat'],'cellAirspace','-v7.3');

%Data with only class B airspace
airspaceB = airspace;
airspaceB.A = zeros(size(airspaceB.A));
airspaceB.B = ones(size(airspaceB.B));
airspaceB.C = zeros(size(airspaceB.C));
airspaceB.D = zeros(size(airspaceB.D));
airspaceB.O = zeros(size(airspaceB.O));
cellAirspace = airspaceB;
save([saveDirectory,'testAirspaceB.mat'],'cellAirspace','-v7.3');

%Data with only class C airspace
airspaceC = airspace;
airspaceC.A = zeros(size(airspaceC.A));
airspaceC.B = zeros(size(airspaceC.B));
airspaceC.C = ones(size(airspaceC.C));
airspaceC.D = zeros(size(airspaceC.D));
airspaceC.O = zeros(size(airspaceC.O));
cellAirspace = airspaceC;
save([saveDirectory,'testAirspaceC.mat'],'cellAirspace','-v7.3');

%Data with only class D airspace
airspaceD = airspace;
airspaceD.A = zeros(size(airspaceD.A));
airspaceD.B = zeros(size(airspaceD.B));
airspaceD.C = zeros(size(airspaceD.C));
airspaceD.D = ones(size(airspaceD.D));
airspaceD.O = zeros(size(airspaceD.O));
cellAirspace = airspaceD;
save([saveDirectory,'testAirspaceD.mat'],'cellAirspace','-v7.3');

%Data with only class O ("other") airspace - i.e., classes E and G
airspaceO = airspace;
airspaceO.A = zeros(size(airspaceO.A));
airspaceO.B = zeros(size(airspaceO.B));
airspaceO.C = zeros(size(airspaceO.C));
airspaceO.D = zeros(size(airspaceO.D));
airspaceO.O = ones(size(airspaceO.O));
cellAirspace = airspaceO;
save([saveDirectory,'testAirspaceO.mat'],'cellAirspace','-v7.3');

%% Temporal Test Data
% Month - TotalTime/MaxOcc of each cell is equal to the corresponding month
[~,M,~,~,~,~] = datevec(double(origCell.Date));
testDataMonth = testData;
testDataMonth.TotalTime = ones(size(origCell.TotalTime));
for m = 1:12
    testDataMonth.TotalTime(ismember(M,m)) = uint32(m);
end
testDataMonth.MaxOcc = testDataMonth.TotalTime;
cell = testDataMonth; 
save([saveDirectory,'testDataMonth.mat'],'cell','-v7.3');
clear testDataMonth; clear M;
clear cell;

% Day - TotalTime/MaxOcc of each cell is equal to the corresponding day of
% week
wd = weekday(origCell.Date);
testDataDay = testData;
testDataDay.TotalTime = ones(size(origCell.TotalTime));
for d = 1:7
    testDataDay.TotalTime(ismember(wd,d)) = uint32(d);
end
testDataDay.MaxOcc = testDataDay.TotalTime;
cell = testDataDay; 
save([saveDirectory,'testDataDay.mat'],'cell','-v7.3');
clear testDataDay; clear wd;
clear cell;

% Hour - TotalTime/MaxOcc of each cell is equal to the corresponding hour
% bin
testDataHour = testData;
for h = 1:8
    testDataHour.TotalTime(ismember(origCell.Time,uint8(h))) = uint32(h);
end
testDataHour.MaxOcc = testDataHour.TotalTime;
cell = testDataHour; 
save([saveDirectory,'testDataHour.mat'],'cell','-v7.3');
clear testDataHour
clear cell;

%% Area Test Data
% Data is only positive for the area of interest
testDataArea = testData;
job = readtable('job_test.dat','ReadVariableNames',false,'Format','%s%s');
job.Properties.VariableNames = {'Name','Value'};
for rr = 1:size(job,1)
    fname = table2cell(job(rr,1));
    fvalue = table2cell(job(rr,2));
    jobdata.(fname{1}) = str2num(fvalue{1}); %#ok<ST2NM>
end

latlim = [48.5,49.5];
lonlim = [-126.5, -125.5];

termaplegend = [ jobdata.BINS_PER_DEGREE jobdata.NORTH_LAT jobdata.WEST_LON];
cellLatCutpoints = fliplr(termaplegend(2):-1/termaplegend(1):termaplegend(2)-jobdata.GRID_Y_NUM/termaplegend(1));
cellLonCutpoints = termaplegend(3):1/termaplegend(1):termaplegend(3)+jobdata.GRID_X_NUM/termaplegend(1);  

% Get bin limits of area/track
[~,~,cellLatLim] = histcounts(latlim,cellLatCutpoints);
[~,~,cellLonLim] = histcounts(lonlim,cellLonCutpoints);

% Get indices into original cell
indmat = reshape(1:jobdata.GRID_Y_NUM*jobdata.GRID_X_NUM,jobdata.GRID_Y_NUM,jobdata.GRID_X_NUM);
geoinds = indmat(cellLatLim(1):cellLatLim(2),cellLonLim(1):cellLonLim(2));
geoinds = geoinds(:);
georeduceinds = ismember(origCell.Cell,uint32(geoinds));
testDataArea.TotalTime(~georeduceinds) = 0;

% Set non-zero data
nonzeroData = uint32([ 1 1 1 1 1 1 1,...                                     
                       1 2 2 2 2 2 1,...
                       1 2 3 3 3 2 1,...
                       1 2 3 4 3 2 1,...
                       1 2 3 3 3 2 1,...
                       1 2 2 2 2 2 1,...
                       1 1 1 1 1 1 1 ]);

for i = 1:numel(geoinds)
    testDataArea.TotalTime(ismember(origCell.Cell,uint32(geoinds(i)))) = nonzeroData(i);
end
                 
% Set MaxOcc and save data
testDataArea.MaxOcc = testDataArea.TotalTime;
cell = testDataArea; 
save([saveDirectory,'testDataArea.mat'],'cell','-v7.3');
clear testDataArea
clear cell; clear origCell;

%% Create Cell Coverage
%Create full cell coverage (100% coverage in each cell)
cellcoverageOrig = ones(10,10,3);
cellcoverage = cellcoverageOrig;
save([saveDirectory,'testCellCoverage.mat'],'cellcoverage','-v7.3');

%Create empty cell coverage (0% coverage in each cell)
cellcoverage = zeros(size(cellcoverageOrig));
save([saveDirectory,'emptyCellCoverage.mat'],'cellcoverage','-v7.3');

%Create half empty cell coverage (100% coverage in half of the cells; 0% coverage in the other half)
cellcoverage = zeros(size(cellcoverageOrig));
cellcoverage(1:2:numel(cellcoverage)) = 1;
save([saveDirectory,'halfEmptyCellCoverage.mat'],'cellcoverage','-v7.3');

%% Cell coverage Test Data
%Create an empty dataset
testDataEmpty = testData;
testDataEmpty.TotalTime = NaN(size(testData.TotalTime));
testDataEmpty.MaxOcc = testDataEmpty.TotalTime;
cell = testDataEmpty;
save([saveDirectory,'emptyTest.mat'],'cell','-v7.3');
clear cell; clear testData;

%Create a half-empty dataset
testDataHalfEmpty = testDataEmpty; clear testDataEmpty;
testDataHalfEmpty.TotalTime(1:2:numel(testDataHalfEmpty.TotalTime)) = uint32(1); 
testDataHalfEmpty.MaxOcc = testDataHalfEmpty.TotalTime;
cell = testDataHalfEmpty;
save([saveDirectory,'halfEmptyTest.mat'],'cell','-v7.3');
clear testDataHalfEmpty

end