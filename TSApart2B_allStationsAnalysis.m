% Time Series Analysis, Part 2B:
% This script reads in station coordinates for ALL stations, 
% calculates trends based on the specified parameters to be estimated in a
% Iteratively Reweighted Least Squares (IRLS) Algorithm
% (polynome degree, oscillations, jumps, ...) and creates output files 
% (plot images, csvs) as a result.
%
%   Input data for this script comprises:
%       - Station Data created by TSA_ReadAndTransform
%           Naming Pattern: <StationName>.mat
%           Content: Table containing 4 columns named "t", "E", "N", "U" and associated data
%       - Jump Database ("Jump Table"):
%           ".csv" file (";" delimiter) with 7 columns named
%           Station = Station Code (e.g. "CONZ", "AREQ", ...)
%           Date = Occurence of Jump in 'yyyy-MM-dd' format
%           Comment = any string
%           Earthquake = 0 (no earthquake) or 1 (is earthquake)
%           HW/SW_Change = 0 (no hardware/software change) or 1 (is hardware/software change)
%           Unknown = 0 (no unknown cause) or 1 (unknown cause)
%           Use = 0 (dont use this jump) or 1 (use this jump)
%
%   Disclaimer: Certain parameters for the (Iteratively Reweighted) Least
%   Squares can be adjusted in the corresponding function file
%   ("computeTrendIRLS.mat")
%
% 
% David Wallinger, DGFI, 2.9.2019

clear variables
close all
addpath('myfunctions')
tic % measure execution time

%% SETTINGS (adapt if necessary) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataStorageLocation = 'station_data'; % Where Station Data (TSA_ReadAndTransform) is stored as ".mat"
stationPositionCSVLocation = 'station_data/_allStationsPosition.csv'; % Simple CSV containing Station Names, Lon, Lat
jumpCSVLocation = 'jumps_version3.csv'; % Location of Jump Table
% jumpCSVLocation = ''; % Location of Jump Table

%%% Trend Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial Trend: Degree (max. 2)
polynDeg = 1;
% polynDeg = 2;

% periods / oscillations in YEARS (=365.25 days), max. 3
P(1) = 1;
P(2) = 1/2;
% P(3) = 6.5;

% Parameter T in [years] for computation of logarithmic transient for
% earthquake events (jumps)
T = 1;

% Model ITRF jumps (set to "true") or ignore ITRF jumps (set to "false")
doITRFjump = false;

%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table (csv) containing ALL results, computed parameters, station
% position, ...
CSVBigTableFilename = 'IRLS_BigTable_AllStations.csv';
% Table (csv) containing only station position and root mean square error
CSVSmallTableFilename = 'IRLS_SmallTable_AllStations.csv';

% Time Series & Trend Plot Image Storage Directory
TSADir = 'stationTSA_results';
imgDir = 'stationTSA_plots';

%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(TSADir, 'dir')
    mkdir(TSADir) % Create Folder if it doesnt exist
    fprintf('TSA Results Storage directory created.\n')
end
% PreProcessing
% load position csv data
stationDataPattern = '%s %f %f %f';
fIDstationData = fopen(stationPositionCSVLocation, 'r');
% get data (station name, lon, lat, h)
stationData = textscan(fIDstationData, stationDataPattern, 'Delimiter',',');
fclose(fIDstationData);

% get station names and cut suffix
fList = dir([dataStorageLocation, '/*.mat']);

% % Debugging: Specific Station
% fList = dir([dataStorageLocation, '/CBSB.mat']);

stationnames = strrep({fList.name}, '.mat', ''); % equals to array of strings of station names
% number of stations (station names) found
nStations = size(stationnames, 2);

% Define Big Table with named columns (*3: E, N and U for every station) % nStations
colNames = {'station','lon','lat','coord','t0','ti','poly0','poly1','poly2', ...
    'osc_T1','oscA1','oscC1','oscS1', ...
    'osc_T2','oscA2','oscC2','oscS2', ...
    'osc_T3','oscA3','oscC3','oscS3', ...
    't_heaviside1','heaviside1','t_heaviside2','heaviside2', ...
    't_heaviside3','heaviside3','t_heaviside4','heaviside4', ...
    't_heaviside5','heaviside5','t_heaviside6','heaviside6', ...
    't_heaviside7','heaviside7','t_heaviside8','heaviside8', ...
    't_heaviside9','heaviside9','t_heaviside10','heaviside10', ...
    't_heaviside11','heaviside11','t_heaviside12','heaviside12', ...
    't_logtrans1','logtrans1', ...
    't_logtrans2','logtrans2', ...
    't_logtrans3','logtrans3', ...
    't_logtrans4','logtrans4', ...
    'rmse'};

CSVBigTable = cell2table(cell(nStations * 3, length(colNames)), 'VariableNames', ...
        colNames);

% Define Small Table with named columns for just rmse
CSVSmallTable = cell2table(cell(nStations, 6), 'VariableNames', ...
    {'station', 'lon', 'lat', 'rmse_E', 'rmse_N', 'rmse_U'});    

%% Loop all found stations
CSVidxCNT = 1; % Index Counter for Big Table

for i = 1:nStations
    
    stationname = stationnames{i};
    fprintf('%s-%d\n',stationname,i)
    % Index for Big Table
    CSVidx = CSVidxCNT:CSVidxCNT + 2; % always 3 rows: E, N and U
    CSVidxCNT = CSVidxCNT + 3;
    
    % Open figure;
    
    %% Call Trend Computation Function for 1 Station
    % Return x struct containing parameters
    x = fTSA_TrendComputation(stationname,jumpCSVLocation, dataStorageLocation, ...
        polynDeg, P, T, doITRFjump, imgDir);  
    
    % Create (0) vectors for parameters
    polynCoeff = zeros(3, 3);
    [osc_T, oscA, oscC, oscS] = deal(zeros(3, 3));
    t_heaviside = strings(1, 12);
    heaviside = zeros(3, 12);
    t_logtrans = strings(1, 4);
    logtrans = zeros(3, 4);
    
    %% fill up parameter vectors with results (metric units: mm!)
    % get metainformation
    t0 = datestr(x.t0, 'yyyy-mm-dd HH:MM'); % first observation
    ti = datestr(x.ti, 'yyyy-mm-dd HH:MM'); % last observation
    
    % get coefficients
    polynCoeff(:, 1:polynDeg + 1) = x.p;
    % get oscillations
    osc_T(:, 1:size(P, 2)) = x.osc_T; % Time (Period in Y = 365.25d)
    oscA(:, 1:size(P, 2)) = x.A; % Amplitude
    oscC(:, 1:size(P, 2)) = x.C; % Cos
    oscS(:, 1:size(P, 2)) = x.S; % Sin
    
    % get heaviside jumps and their time of occurence
    % events
    t_heaviside(1:size(x.t_heaviside(1, :), 2)) = ...
        datestr(x.t_heaviside(1, :), 'yyyy-mm-dd HH:MM');
    
    % jumps
    heaviside(:, 1:size(x(1).heaviside, 2)) = x.heaviside;
    
    % get logarithmic transients and their time of occurence (->EQ event)
    % events
    t_logtrans(1:size(x.t_logtrans(1, :), 2)) = ...
        datestr(x.t_logtrans(1, :), 'yyyy-mm-dd HH:MM');
    
    % transients
    logtrans(:, 1:size(x.t_logtrans(1, :), 2)) = x.logtrans;
    
    %% append to big table
    
    % Stationname
    CSVStationnameCol = cell(3, 1);
    CSVStationnameCol(:) = {stationname};
    CSVBigTable(CSVidx, 'station') = CSVStationnameCol;
    % lon, lat
    CSVBigTable{CSVidx, 'lon'} = num2cell(repmat(stationData{3}(i), 3, 1));
    CSVBigTable{CSVidx, 'lat'} = num2cell(repmat(stationData{2}(i), 3, 1));
    
    % t0,ti (needs conversion to cell)
    CSVBigTable(CSVidx, 't0') = cellstr(t0);
    CSVBigTable(CSVidx, 'ti') = cellstr(ti);
    
    % coord
    CSVBigTable(CSVidx, 'coord') = cellstr(x.coordinate);
    
    % polynome coefficients
    CSVBigTable{CSVidx, {'poly0', 'poly1', 'poly2'}} = num2cell(polynCoeff);
        
    % oscillation event datetimes
    CSVBigTable{CSVidx, {'osc_T1', 'osc_T2', 'osc_T3'}} = num2cell(osc_T);
    % oscillation parameters
    CSVBigTable{CSVidx, {'oscA1', 'oscA2', 'oscA3'}} = num2cell(oscA);
    CSVBigTable{CSVidx, {'oscC1', 'oscC2', 'oscC3'}} = num2cell(oscC);
    CSVBigTable{CSVidx, {'oscS1', 'oscS2', 'oscS3'}} = num2cell(oscS);
    
    % heaviside event datetimes (columns 22:2:36), same datetime x3
    if length(t_heaviside) <= 12
        CSVBigTable{CSVidx, {...
            't_heaviside1', 't_heaviside2', 't_heaviside3', ...
            't_heaviside4', 't_heaviside5', 't_heaviside6', ...
            't_heaviside7', 't_heaviside8', 't_heaviside9', ...
            't_heaviside10', 't_heaviside11', 't_heaviside12'}} = ...
            cellstr([t_heaviside; t_heaviside; t_heaviside]);
        
        % heaviside jumps (columns 23:2:37)
        CSVBigTable{CSVidx, {...
            'heaviside1', 'heaviside2', 'heaviside3', ...
            'heaviside4', 'heaviside5', 'heaviside6', ...
            'heaviside7', 'heaviside8', 'heaviside9', ...
            'heaviside10', 'heaviside11', 'heaviside12'}} = num2cell(heaviside);
    else
        fprintf('Not all Jumps stored, because total number of jumps exceeded 12.\n')
        % Not more than 12 Jumps can be stored, from 13 onwards will not be
        % stored!
        CSVBigTable{CSVidx, {...
            't_heaviside1', 't_heaviside2', 't_heaviside3', ...
            't_heaviside4', 't_heaviside5', 't_heaviside6', ...
            't_heaviside7', 't_heaviside8', 't_heaviside9', ...
            't_heaviside10', 't_heaviside11', 't_heaviside12'}} = ...
            cellstr([t_heaviside(1:12); t_heaviside(1:12); t_heaviside(1:12)]);
        
        % heaviside jumps (columns 23:2:37)
        CSVBigTable{CSVidx, {...
            'heaviside1', 'heaviside2', 'heaviside3', ...
            'heaviside4', 'heaviside5', 'heaviside6', ...
            'heaviside7', 'heaviside8', 'heaviside9', ...
            'heaviside10', 'heaviside11', 'heaviside12'}} = num2cell(heaviside(:, 1:12));
        
    end
    % eq event datetimes, same datetime x3
    CSVBigTable{CSVidx, {'t_logtrans1', 't_logtrans2', 't_logtrans3', 't_logtrans4'}} = cellstr([t_logtrans; t_logtrans; t_logtrans]);
    
    % eq logarithmic transients
    CSVBigTable{CSVidx, {'logtrans1', 'logtrans2', 'logtrans3', 'logtrans4'}} = num2cell(logtrans);
    
    % rsme
    CSVBigTable{CSVidx, 'rmse'} = num2cell(x.rmse);
    
    %% append to small table
    CSVSmallTable{i, 'station'} = {stationname};
    CSVSmallTable{i, 'lon'} = {stationData{3}(i)};
    CSVSmallTable{i, 'lat'} = {stationData{2}(i)};
    CSVSmallTable{i, 'rmse_E'} = {x.rmse(1)}; % 1st element: E
    CSVSmallTable{i, 'rmse_N'} = {x.rmse(2)}; % 1st element: N
    CSVSmallTable{i, 'rmse_U'} = {x.rmse(3)}; % 1st element: U
    
end

% write table to csv file
writetable(CSVBigTable,fullfile(TSADir, CSVBigTableFilename),'writevariablenames',1)
writetable(CSVSmallTable,fullfile(TSADir, CSVSmallTableFilename),'writevariablenames',1)

toc % Print Runtime

fprintf('Finished!\n');
