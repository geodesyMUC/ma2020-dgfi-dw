clear variables; close all;

% Visualizes different trends of the same station overlying their original timeline
% Input is the result csv file computed by TSApart2A_TrendComputation
%

addpath('myfunctions')

%% INPUT data settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputFolder = 'station_data'; % Where Station Data (TSA_ReadAndTransform) is stored as ".mat"
jumpCSVLocation = 'jumps_version3.csv'; % Location of Jump Table/Jump Database
resultsLocation = 'TSA_TrendComputationResults';

%%% Name of station to be analysed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTION FOR ANALYSIS

% stationname = 'MEXI'; % Mexicali, Mexico
% stationname = 'ALAR'; % Arapiraca, Brazil
stationname = 'AREQ'; % Arequipa, Peru
% stationname = 'CONZ'; % Concepcion, Chile
% stationname = 'OAX2'; % Oaxaca, Mexico
% stationname = 'CUEC'; % Cuenca, Ecuador
% stationname = 'MZAE'; % Santa Rosa, Argentina (missing jump)
% stationname = 'NEIL'; % Ciudad Neilly, Costa Rica
% stationname = 'RWSN'; % Rawson, Argentina

% result files to be analysed
% get files according to their filename in the specified directory
[Results, Result_params, Results_names] = findStationData(stationname, resultsLocation, 3);

% other variables (constants)
coordinateSTR = {'East', 'North', 'Up'}; % used to label plots
titleString = ['Comparison of LSE Results for "', stationname, '"'];

%% Add ITRF Realization Change Jumps (for plots)
% format yyyy-MM-dd,
jumpITRF(1) = datetime('2011-04-17', 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
jumpITRF(2) = datetime('2017-01-29', 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
% can be extended! keep format "yyyy-MM-dd", time will be 0:00am)

%% Select Station Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look in matched files
% naming pattern for station measurement data has to be "<StationName>.mat"
fpathpattern = fullfile(inputFolder, sprintf('%s.mat', stationname)); % match pattern
if exist(fpathpattern, 'file')
    load(fpathpattern)
    fprintf('Evaluating Station "%s".\n', stationname);
else
    error('Specified station "%s" - associated .mat-file could not be found!', stationname)
end

% Reassign currStation variable
data = currStation.Data;
dataSTATION_NAME = currStation.Station;

% Load Jump Table using custom import function (read from database)
dataJump = importfileJumpCSV(jumpCSVLocation);

% Get jumps for this station
TFJump = strcmp(dataJump{:, 'Station'}, stationname); % selection logical
fprintf('%d jumps for this station found. (excluding ITRF changes)\n', nnz(TFJump));

% All Jumps for current station
currStationJumps = dataJump(TFJump, :);
% All Jumps for current station where "Use" Column Element is set to "1"
jumpTable = currStationJumps(currStationJumps.Use == 1, :);

% Adjust plot names
result_names_cell = cellfun(@(x) strrep(x, stationname, ''), Results_names, 'UniformOutput', false);
result_names_cell = cellfun(@(x) strrep(x, '_', ' '), result_names_cell, 'UniformOutput', false);
result_names_cell = cellfun(@(x) strrep(x, '.csv', ''), result_names_cell, 'UniformOutput', false);

%% Visualize E - N - U time series separately
for i = 1:3
    
    figure
    
    % Plot Measurements
    pPts = plot(data{:, 'date'}, data{:, 2 + i}, ...
        '.', 'markersize', 4, ...
        'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Observations');
    
    ax = gca;
    y1 = ax.YLim(1); % axis MIN
    y2 = ax.YLim(2); % axis MAX
    
    grid on
    hold on
    
    % Plot Trend (j = 1:number of files)
    for j = 1:size(Results, 1)
%         % set up plot name
%         pTrendSplitName = strrep(Results_names{j}, '_', '-');
%         pTrendSplitName = strsplit(pTrendSplitName, '.');
%         pTrendName = [pTrendSplitName{1}, '.', pTrendSplitName{2}];
        % plot
        pTrend = plot(Results{j, 1}, Results{j, 1 + i}, 'DisplayName', result_names_cell{j});
    end
    
    ylabel(sprintf('%s [mm]', coordinateSTR{i}))
    xlim([min(data{:, 'date'}) max(data{:, 'date'})])
    title(titleString);
    
    % Plot Jumps from Jump Table
    for j = 1:size(jumpTable, 1)
        if jumpTable{j, 4} > 0 && jumpTable{j, 2} > data{1, 'date'} % if type matches and jump is AFTER first obs.
            % Earthquake
            pEq = plot([jumpTable{j, 2}; jumpTable{j, 2}], [y1; y2], ...
                'color', [102, 51, 0]./255, 'DisplayName', 'EQ Jump');
        elseif jumpTable{j, 5} > 0 && jumpTable{j, 2} > data{1, 'date'}
            % Antenna Change, ...
            pAnt = plot([jumpTable{j, 2}; jumpTable{j, 2}], [y1; y2], ...
                'color', [0.75, 0, 0.75], 'DisplayName', 'HW Change');
        elseif jumpTable{j, 6} > 0  && jumpTable{j, 2} > data{1, 'date'}
            % Unknown cause
            pUnk = plot([jumpTable{j, 2}; jumpTable{j, 2}], [y1; y2], ...
                'color', [0.4660, 0.6740, 0.1880], 'DisplayName', 'Unkn.Jump');
        end
    end
    
    % plot itrf jump vertical lines
    for j = 1:length(jumpITRF)
        pITRF = plot([jumpITRF(j); jumpITRF(j)], [y1; y2], '--', ...
            'color', [160, 160, 200]./255, 'DisplayName', 'new ITRF');
    end
    
    hold off
    legend show
    % set location of legend to "best", needs legend handle
    hLegend = findobj(gcf, 'Type', 'Legend');
    hLegend.Location = 'northeastoutside';
    % change figure dimensions
    set(gcf, 'InnerPosition', [1000 1000 1500 600]); % large figure
end

%% Visualize WRMS/RMS as bar graph
fprintf('Evaluating Parameters of "%s".\n', stationname);

for i = 1:3 % E, N, U
    c = categorical(result_names_cell);
    result_param_array = cell2mat(Result_params);
    plotOffsetLimit = 0.1;
    
    figure;
    subplot(2, 1, 1) % RMS bar plot
    data = result_param_array(2:3:end, i + 1);
    [data, I] = sort(data, 'descend');
    c = categorical(result_names_cell(I), result_names_cell(I));
    barh(c, data, 'FaceColor', [0 0.4470 0.7410]);
    grid on
    title(sprintf('Station "%s", RMS Comparison(%s)', stationname, coordinateSTR{i}));
    ylabel('error [mm]')
    xlim([min(data)-plotOffsetLimit, max(data)+plotOffsetLimit])
    
    subplot(2, 1, 2) % WRMS bar plot
    data = result_param_array(3:3:end, i + 1);
    [data, I] = sort(data, 'descend');
    c = categorical(result_names_cell(I), result_names_cell(I));  
    barh(c, data, 'FaceColor', [0.8500 0.3250 0.0980]);
    grid on
    title(sprintf('Station "%s", WRMS Comparison(%s)', stationname, coordinateSTR{i}));
    ylabel('error [mm]')
    xlim([min(data)-plotOffsetLimit, max(data)+plotOffsetLimit] )
    
end

function [result_array, result_params, result_names] = findStationData(stationname, location, headerCnt)
% get file list
fList = dir([location, '/', stationname, '*.csv']);
% fList([1 5 6 7]) = [] % REMOVE FILES FROM ANALYSIS MANUALLY
fprintf('Files:\n%s', sprintf('%s\n', fList.name));
n_files = length(fList);
% get data from result files (result files (datetime unix, E, N, U)
% preallocating
result_array = cell(n_files, 4);
result_names = cell(n_files, 1);
result_params = cell(n_files, 1); 

for i = 1:n_files
    readInData = csvread(fullfile(fList(i).folder, fList(i).name));
    
    result_params{i} = readInData(1:headerCnt, :); % first 3 rows
    
    % remaining rows -> data
    result_array{i, 1} = datetime(readInData(headerCnt + 1:end, 1), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
    result_array{i, 2} = readInData(headerCnt + 1:end, 2); % E
    result_array{i, 3} = readInData(headerCnt + 1:end, 3); % N
    result_array{i, 4} = readInData(headerCnt + 1:end, 4); % U
    
    % filename
    result_names{i} = fList(i).name;
end

end