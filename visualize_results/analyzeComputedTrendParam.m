% This script creates plots of the specified parameter (e.g. RMS)
% Input is the BIG table computed by "TSApart2B"
%
%   The idea is to get quick insights on how the resulting station RMS
%   errors are distributed, so e.g. outliers can be identified.
%
% David Wallinger, DGFI, August 2019

clear variables
close all
addpath('myfunctions') % Add Function Storage to PATH

%% SETTINGS (adapt if necessary) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BigTableLocation = 'stationTSA_results/IRLS_BigTable_AllStations.csv';

%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = importStationTrendDataCSV(BigTableLocation);
% data = importAllStationTrendData('TrendParam_AllStations.csv');

% Split all Stations East North Up
dataE = data(strcmp(data.coord, 'E'), :);
dataN = data(strcmp(data.coord, 'N'), :);
dataU = data(strcmp(data.coord, 'U'), :);

%% visualization #1
figure
title('Root Mean Square Error of Trend for ENU')
hold all
plot(dataE.rmse, ones(length(dataE.rmse), 1) .* 1, 'r+')
plot(dataN.rmse, ones(length(dataN.rmse), 1) .* 2, 'b+')
plot(dataU.rmse, ones(length(dataU.rmse), 1) .* 3, 'g+')

% text labels
text(dataE.rmse, ones(length(dataE.rmse), 1) .* 1 + 0.08, dataE.station, ...
    'FontSize',7, 'Rotation', 45);
text(dataN.rmse, ones(length(dataN.rmse), 1) .* 2 + 0.08, dataN.station, ...
    'FontSize',7, 'Rotation', 45);
text(dataU.rmse, ones(length(dataU.rmse), 1) .* 3 + 0.08, dataU.station, ...
    'FontSize',7, 'Rotation', 45);

ylim([0 4])
yticks([1,2,3])
yticklabels({'E', 'N', 'U'})
xlabel('RMSE [mm]')
grid on
hold off

%% visualization #2 - parameter evaluation
figure
title('Box Plot of RMSE of computed trends')
boxplot([dataE.rmse dataN.rmse dataU.rmse], ...
    'Labels',{'East','North', 'Up'})
ylabel('Root Mean Square Error [mm]')

