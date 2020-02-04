% Time Series Analysis, Part 1:
% This script reads in a ".neh" file containing station observations and creates
% Plots and Data for further processing and evaluation
%
%   The Input file must contain:
%   DATE STATION NAME PHI LAMBDA ELEVATION (White Space Delimiter)
%
%   The data is being converted into North, East, Up. A ".mat" file for
%   every station containing the converted coordinates will be created.
%
%   A folder "temp" will be created, containing a ".mat" file with all
%   station positions for faster processing when running this script
%   multiple times. THIS FILE NEEDS TO BE DELETED IF NEW DATA IS TO BE
%   PROCESSED!!!
%
%   David Wallinger, DGFI, 1.7.2019

clear variables
close all
addpath('myfunctions') % Add Function Storage to PATH
%% SETTINGS (adapt if necessary) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input File Location & Name 
fname = 'raw_data/SirgasCoord.neh';

% Output Settings (changing them is not necessary)
outputFolder = 'station_data'; % Storage for station coordinate files
outputFilename1 = '.mat'; % File Extension (DO NOT CHANGE, HAS TO BE ".mat")
outputFilename2 = '-TSvis.png'; % Station TS Plot Image Naming Pattern. 
% -> Station name will be attached, e.g. "CONZ-TSvis.png"
outputStationPositionFile = '_allStationsPosition.csv'; % Name of simple csv file
% containing: Station Name, Station Latitude, Station Longitude and Station Elevation

%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create output folder for data and images
if ~exist(outputFolder, 'dir')
   mkdir(outputFolder) 
end
% Create temp folder for script data storage
if ~exist('temp', 'dir')
  mkdir('temp');
end

fID = fopen(fullfile(outputFolder, outputStationPositionFile), 'wt');

% other variables
coordinateSTR = {'E', 'N', 'U'};

% Check if file already exists. (Saves Time when running again)
% If not, use read in function
if ~exist(fullfile('temp', 'dataAllStations.mat'), 'file')
    % read in
    data = importfileBLh(fname);
    % save
    save(fullfile('temp', 'dataAllStations.mat'), 'data');
else
    % If exists -> just load .mat file
    load(fullfile('temp', 'dataAllStations.mat'), 'data');
    fprintf('Using existing file for station data!\n');
end

% create txts and png TS plots
nStations = length(data);

for i = 1:nStations
    % current station information
    currStation.Station = data(i).Station{:};
    currStation.Data = data(i).Data;
    
    filepath1 = fullfile(outputFolder, [currStation.Station, outputFilename1]);
    filepath2 = fullfile(outputFolder, [currStation.Station, outputFilename2]);
    
    % Save ".mat" file
    save(filepath1, 'currStation');
    
    % Create TS Print (image)
    % VisualizeTS_ENU function needed!
    figTS = figure('visible','off');
    VisualizeTS_ENU(currStation.Data, currStation.Station);
    
    % Save image file
    saveas(figTS, filepath2);
    close(figTS)
    
    % Print Station Position to ".csv"
    fprintf(fID, '%s,%f,%f,%f\n', currStation.Station, ...
        data(i).StationPosition1(1), data(i).StationPosition1(2), data(i).StationPosition1(3));
end

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OUTPUT = importfileBLh(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   SIRGASCOORD = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   SIRGASCOORD = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   SirgasCoord = importfile('SirgasCoord.neh', 5, 165241);
%
%    See also TEXTSCAN.

%% Initialize variables.
if nargin<=2
    startRow = 5;
    endRow = inf;
end

%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%10s%5s%19f%15f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Remove white space around all cell columns.
dataArray{1} = strtrim(dataArray{1});
dataArray{2} = strtrim(dataArray{2});

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.
%% STATIONS
% Get stations
[stations, stationsIDX1, stationsIDX] = unique(dataArray{2});
nStations = length(stations);
% % Count measurements per station
% stationsCNT = accumarray(stationsIDX, 1);

% Use all stations
stationIDX = 1:nStations;

%% Get data for all stations
% parameters
spheroid = referenceEllipsoid('GRS 80');
strvec120000= ' 12:00:00'; % HH:mm:ss
% preallocate Output
OUTPUT = struct('Station', cell(1, nStations), 'Data', cell(1, nStations));

for i = 1:nStations
    
    TF = stationsIDX == stationIDX(i); % current station logical
    sIDX1 = stationsIDX1(i); % current station first idx in vec
    
    %% Transform
    % [xEast,yNorth,zUp] = geodetic2enu(lat,lon,h,lat0,lon0,h0,spheroid)
    % get origin (use 1!)
    
    [E, N, U] = geodetic2enu(... 
        dataArray{3}(TF),... % Latitude
        dataArray{4}(TF),... % Longitude
        dataArray{5}(TF), ... % Ell. Height
        dataArray{3}(sIDX1), ... % Origin Latitude (use 1. obs)
        dataArray{4}(sIDX1), ... % Origin Longitude (use 1. obs)
        dataArray{5}(sIDX1), ... % Origin Height (use 1. obs)
        spheroid, 'degrees'); % other parameters
       
    %% Station datetime    
    % add 12:00 string to adjust time!
    % datetime vector
    D = datetime(strcat(dataArray{1}(TF), strvec120000), ...
        'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone','UTC');
    reftime = D(1); % first date as a reference
    % create vector for seconds (intervals) -> compute difference with
    % etime
    DT_s = etime(datevec(D), repmat(datevec(reftime), ...
        size(D, 1)));
    
    %% Output
    % Reassign & set up new data table
%     dataSTATION_NAME = stations(stationIDX(i));
    OUTPUT(i).Data = table(...
        D, ... 
        DT_s, ... % second intervals
        E .* 1e3, ... % mm!
        N .* 1e3, ... % mm!
        U .* 1e3, ... % mm!
        'VariableNames', ...
        {'date', 't', 'E', 'N', 'U'}...
        );
    OUTPUT(i).Station = stations(stationIDX(i));
    % Phi Lambda h
    OUTPUT(i).StationPosition1 = [...
        dataArray{3}(sIDX1), dataArray{4}(sIDX1), dataArray{5}(sIDX1)];
end
end
