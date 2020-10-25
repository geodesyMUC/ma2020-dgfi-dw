% This script creates map plots for the examined region 
% (determined by station coordinates).
% Input is the BIG table computed by "TSApart2B"
% After specifying the BIG Table CSV Location :
% 
% First, a certain parameter from the BIG Table csv needs to be specified.
% E.g.:
%   plotDataC = dataU.oscA1; (Marker Color)
%   plotDataS = 20; (Marker Size, can be either a parameter or static)
%       -> Select dataE dataN or dataU
% Secondly, a certain date (yyyy-mm-dd) can be specified. The script uses the Tectonic
% Plate Boundaries from the corresponding directory (needs to be specified
% too) to create a map scheme. Then, the examined jumps will being visualized
% depending on the station positions. 
% 
% Use the data cursor to get
% information about the selected station!
% 
% Input is the BIG table computed by "TSApart2B"
% David Wallinger, DGFI, August 2019

clear variables
close all
addpath('myfunctions') % Add Function Storage to PATH

%% SETTINGS (adapt if necessary) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BigTableLocation = 'stationTSA_results/IRLS_BigTable_AllStations.csv'; % Big Table Location

plateBoundariesDir = 'TectonicPlateBndsLatLon'; % Plate Boundaries Location

%%% Do not change !!!
data = importStationTrendDataCSV(BigTableLocation); % Import Big Table
% Split all Stations East North Up
dataE = data(strcmp(data.coord, 'E'), :);
dataN = data(strcmp(data.coord, 'N'), :);
dataU = data(strcmp(data.coord, 'U'), :);
%%% !!! %%%

%%%% Select data for evaluation (COLOR SCALE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotDataC = dataE.poly1;
% plotDataC = dataU.oscA2;
plotDataC = dataU.oscA1;
% plotDataC = sqrt(dataE.rmse.^2 + dataN.rmse.^2);
% plotDataC = dataU.rmse.^2;

%%%% Select data for evaluation (SIZE SCALE, can be scalar (=static)) %%%%%
% plotDataS = 20; % static
% plotDataS = dataE.poly1;
% plotDataS = dataU.oscA2;
plotDataS = dataU.oscA1;
% plotDataS = sqrt(dataE.rmse.^2 + dataN.rmse.^2);
% plotDataS = dataU.rmse.^2;

%%%% Select Jumps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jumpOccurenceDatetime = '2010-02-27'; % Chile EQ
% jumpOccurenceDatetime = '2016-04-16'; % Equador EQ
% jumpOccurenceDatetime = '2015-09-16'; % 
jumpOccurenceDatetime = '2014-04-12';

%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualization #1 - Trends
close all % !!!
% createCustomMap(dataE.lat, dataE.lon)

% Reassign Position of Station
stationLat = dataE.lat;
stationLon = dataE.lon;

% Create WorldMap
figure;

worldmapBuffer = 5; % deg
worldmap(...
    [min(stationLat), max(stationLat) + worldmapBuffer], ...
    [min(stationLon) - worldmapBuffer, max(stationLon) + worldmapBuffer]);
load coastlines
plotm(coastlat, coastlon, 'k')
hold on

% Read Plate Bndaries
rawPlateFileList = dir([plateBoundariesDir '/*']);
rawPlateCoordList = cell(length(rawPlateFileList), 2);

for i = 1:length(rawPlateFileList)
    if regexpi(rawPlateFileList(i).name, '^[a-z]{2}$')
        [plateLat, plateLon] = importTectonic(fullfile(rawPlateFileList(1).folder, rawPlateFileList(i).name));
        
        plotm(plateLat, plateLon, 'color', [180 180 100]./255);
        
        rawPlateCoordList{i , 1} = plateLat;
        rawPlateCoordList{i , 2} = plateLon;
        
    end
end

plateCoordList = rawPlateCoordList(~cellfun('isempty',rawPlateCoordList));
plateCoordList = reshape(plateCoordList, length(plateCoordList)/2, 2);

% Data Plot
plotAllScatter = scatterm(stationLat, stationLon, plotDataS, plotDataC, 'filled', 'o');

colorbar
colormap jet

figure
hold on
plot(plotDataC, ones(length(plotDataC)), '+'); % simple 

% text labels
text(plotDataC, ones(length(stationLat), 1) + 0.05, dataE.station, ...
    'FontSize',7, 'Rotation', 90);

scatter(plotDataC, ones(size(plotDataC, 1), 1) .* 2,  ...
    ones(size(plotDataC, 1), 1) .* plotDataS, plotDataC, '+'); % colored/scaled
hold off
title('Data Distribution')
ylim([0.5 2.5])
yticks([1 2])
grid on
colorbar
colormap jet

%% visualization #2 - jump evaluation

jumpMagnitude = zeros(length(stationLat), 3); % E N U

% Loop and get matches
for i = 22:2:44
    jumpMatchLogical = dataE{:, i} == jumpOccurenceDatetime;
    jumpMatchIndices = find(jumpMatchLogical == 1);
    % Assign E N U
    jumpMagnitude(jumpMatchIndices, 1) = dataE{jumpMatchIndices, i + 1}; % next element correspondends to jump value
    jumpMagnitude(jumpMatchIndices, 2) = dataN{jumpMatchIndices, i + 1}; % next element correspondends to jump value
    jumpMagnitude(jumpMatchIndices, 3) = dataU{jumpMatchIndices, i + 1}; % next element correspondends to jump value   
end


% Loop E N U to get plots
labelCell = {'East', 'North', 'Up'};
for i = 1:3
    
    isJump = abs(jumpMagnitude(:, i)) > 0; % Logical for relevant jumps
    
    figure%('DeleteFcn','doc datacursormode');
    % Custom Datacursor
    dcm_obj = datacursormode(gcf);

    title(['Jumps in ', labelCell{i}, ' component, Jump Event: ', jumpOccurenceDatetime])
    
%     % Custom Datacursor
%     dcm_obj = datacursormode(gcf);
%     set(dcm_obj,'UpdateFcn',{@myupdatefcn, [plotAllPts.Children.XData', plotAllPts.Children.YData'], ...
%         dataE.station, jumpMagnitude(:, i)})

    % set colormap
    if max(jumpMagnitude(:, i)) > abs(min(jumpMagnitude(:, i)))
        colormap(flipud(autumn))
    else
        colormap(autumn)
    end
    colorbar
        
    worldmap(...
        [min(stationLat), max(stationLat) + worldmapBuffer], ...
        [min(stationLon) - worldmapBuffer, max(stationLon) + worldmapBuffer]);
    mlabel('off')
    gridm('off')
    % Plot Plate Boundaries
    for j = 1:length(plateCoordList)
        plotm(plateCoordList{j, 1}, plateCoordList{j, 2}, 'Color', [180 180 100]./255)
    end
    
    % Plot Coastlines
    hold on
    plotm(coastlat, coastlon, 'k')
    plotAllPlot = plotm(stationLat, stationLon, 'o', 'Color', [170 170 170]./255, 'MarkerSize', 3);

    hold on
    % Add the data ---
    scatterm(stationLat(isJump), stationLon(isJump), ...
        35, jumpMagnitude(isJump, i), 'o', 'filled', ...
        'MarkerEdgeColor',[0 0 0])
    
    % Set Datacurser function
    set(dcm_obj,'UpdateFcn',{@myupdatefcn, [plotAllPlot.XData', plotAllPlot.YData'], ...
        [plotAllScatter.Children.XData', plotAllScatter.Children.YData'], ...
        dataE.station, jumpMagnitude(:, i)})
    refresh
    
end
fprintf('Use Data Cursor in figure to get Station Information!\n')