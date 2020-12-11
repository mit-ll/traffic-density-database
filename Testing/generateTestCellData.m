% Copyright 2008 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%
%The purpose of this script is to generate simple test data that can be
%used in the TrafficDensity tool unit tests. 

% One Year of Data
startDate = datenum([2018,1,1]);
endDate = datenum([2018,12,31]);

timeValues = 1:8; % Eight 3-hour bins

cellValues = 1:100; %10x10 test area --> 100 cells

heightBins = 0:2; % 3 altitude bins - 0-FL100, FL100-FL200, FL200-FL300

ACCategory = 0:1; % 0 = discrete-code; 1 = 1200-code. In the test data, 1200-code aircraft in height/altitude bin 0. Discrete-code aircraft in all height/altitude bins

% Initialize data fields
Date = uint32.empty(1168000,0); 
Time = uint8.empty(1168000,0); 
Cell = uint32.empty(1168000,0);
Height = uint32.empty(1168000,0); 
Category = uint32.empty(1168000,0); 
TotalTime = uint32.empty(1168000,0); 
MaxOcc = uint32.empty(1168000,0); 

% Loop through values to populate table for discrete-code aircraft
count = 1;
for date = startDate:endDate % Date
    for t = timeValues % Time
        for c = cellValues % Cell
            for h = heightBins % Height
                Date(count) = date;
                Time(count) = uint8(t);
                Cell(count) = uint32(c);
                Height(count) = uint32(h);
                Category(count) = uint32(0); % Discrete-code
                TotalTime(count) = uint32(1); % Initialize as 1's for now. Will be manipulated in generateTestData.m
                MaxOcc(count) = uint32(1); % Initialize as 1's for now. Will be manipulated in generateTestData.m
                count = count + 1;
            end
        end
    end
end
        
% Loop through values to populate table for 1200-code aircraft
for date = startDate:endDate % Date
    for t = timeValues % Time 
        for c = cellValues % Cell
            for h = 0 %Only populate in the first altitude bin
                Date(count) = date;
                Time(count) = uint8(t);
                Cell(count) = uint32(c);
                Height(count) = uint32(h);
                Category(count) = uint32(1); % 1200-code
                TotalTime(count) = uint32(1); % Initialize as 1's for now. Will be manipulated in generateTestData.m
                MaxOcc(count) = uint32(1); % Initialize as 1's for now. Will be manipulated in generateTestData.m
                count = count + 1;
            end
        end
    end
end

% Create Table
cell = table(Date', Time', Cell', Height', Category', TotalTime', MaxOcc',...
    'VariableNames',{'Date','Time','Cell','Height','Category','TotalTime','MaxOcc'});

% Save Data
save([getenv('TrafficDensityPath') filesep 'Testing' filesep 'TestData' filesep 'cellOrigTest.mat'],'cell');