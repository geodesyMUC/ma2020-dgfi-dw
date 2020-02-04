% Time Series Analysis, Part 2A: MASTER THESIS Nov2019
% This script reads in station coordinates for a SINGLE station, 
% calculates a trend based on the specified parameters to be estimated in a
% Iteratively Reweighted Least Squares (IRLS) Algorithm
% (polynome degree, oscillations, jumps, ...) and plots the result.
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
% David Wallinger, DGFI, 5.8.2019

clear variables;
close all;
addpath('myfunctions')

%% SETTINGS (adapt if necessary) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT data settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputFolder = 'station_data'; % Where Station Data (TSA_ReadAndTransform) is stored as ".mat"
% jumpCSVLocation = 'src/jumps_version3.csv'; % Location of Jump Table/Jump Database
% itrf_changes_textfile = 'src/itrf_changes.txt';

inputFolder = 'station_data_dailyXYZfiles'; % Where Station Data (TSA_ReadAndTransform) is stored as ".mat"
jumpCSVLocation = 'src/jumps_dailyXYZfiles.csv'; % Location of Jump Table/Jump Database
itrf_changes_textfile = 'src/itrf_changes.txt';

%%% Name of station to be analysed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTION FOR ANALYSIS

% stationname = 'MEXI'; % Mexicali, Mexico
% stationname = 'ALAR'; % Arapiraca, Brazil
% stationname = 'AREQ'; % Arequipa, Peru
% stationname = 'CONZ'; % Concepcion, Chile
% stationname = 'OAX2'; % Oaxaca, Mexico
% stationname = 'CUEC'; % Cuenca, Ecuador
% stationname = 'MZAE'; % Santa Rosa, Argentina (missing jump)
% stationname = 'NEIL'; % Ciudad Neilly, Costa Rica
% stationname = 'RWSN'; % Rawson, Argentina
% stationname = 'PBJP';

% stationname = '21701S007A03'; % KSMV
stationname = '21702M002A07'; % MIZU

%%% Trend Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial Trend: Degree
polynDeg = 1; % integer degree number
% polynDeg = 2;
% polynDeg = 3;

% periods / oscillations in YEARS (=365.25 days) in vector form
P = [];
% P(1) = 1;
% P(2) = 1/2;
% P(3) = 10;

% convert oscillations to angular velocity
W = 2 * pi ./ P;

% Parameter T in [years] for computation of logarithmic transient for
% earthquake events (jumps)
T = 1;

% Model ITRF jumps (set to "true") or ignore ITRF jumps (set to "false")
doITRFjump = false; % E - N - U

% Additional Parameters for LSE/IRLSE (can be adjusted with care)
KK = 0; % n of iterations for IRLS
p = 2.0; % L_p Norm for IRLS
outl_factor = 4; % median(error) + standard deviation * factor -> outlier

%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logFileFolder = 'TSA_TrendComputationResults'; % output: log file directory
logFile = [stationname, '_TrendComputation_log.txt']; % output: log file name

%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(logFileFolder, 'dir')
   mkdir(logFileFolder)
   fprintf('Result file storage directory "%s" created.\n', logFileFolder);
end

% Open Log File and get identifier
fID = fopen(fullfile(logFileFolder, logFile), 'wt');

% other variables (constants)
coordinateSTR = {'E', 'N', 'U'}; % used to label plots
jump_category_names = {'Earthquake', 'SW/HW-Change', 'Unknown'}; % corresponds to jump table columns

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

% Load Jump Table using custom import function
dataJump = importfileJumpCSV(jumpCSVLocation);

% Get jumps for this station
TFJump = strcmp(dataJump{:, 'Station'}, stationname); % selection logical
fprintf('%d jumps for this station found.\n', nnz(TFJump));

% All Jumps for current station
currStationJumps = dataJump(TFJump, :);
% All Jumps for current station where "Use" Column Element is set to "1"
currStationJumps = currStationJumps(currStationJumps.Use == 1, :);

%% Visualise Raw Time Series and Jumps from Jump Table
% Plot TS and Jumps
figVis = figure;
% VisualizeTS_ENU(data, dataSTATION_NAME, currStationJumps{:, 2})
VisualizeTS_ENU2(data, dataSTATION_NAME, currStationJumps{:, 2}, ...
    'jumpTypes', currStationJumps{:, 4:5});

%% PREPARE PARAMETERS FOR TREND ESTIMATION
t = data{:, 't'}; % [seconds]; for years, do /(365.25 * 86400);
t0 = data{1, 'date'}; % beginning of ts as datetime

%% Jump table(s)
HJumps = getRelativeJumps(currStationJumps{:, 2}, t0);

% Distinguish between EQs (invokes log. transient) and other jumps (unknown cause or HW
% change -> do not invoke transient)
EQLogical = logical(currStationJumps{:, 4}); % 1:=earthquake; ~1:=no earthquake
% Get Earthquake Jump Vector - Those Jumps will invoke a log. transient
EQJump = getRelativeJumps_eq(currStationJumps{:, 2}, t0, EQLogical); 

%% Add ITRF Realization Change Jumps
if doITRFjump
    fprintf('Considering ITRF jumps, doITRFjump is set to "true".\n')
    
    jumps0itrf = getRelativeITRFJumps(t0, itrf_changes_textfile);
    
    % append itrf jumps to heaviside jumps vector
    HJumps = [HJumps; jumps0itrf];
else
    fprintf('Not considering ITRF jumps, doITRFjump is set to "false".\n')
end

%% Prepare Trend Estimation
% preallocate cell arrays
result_parameters = cell(3, 2);
outlier_logicals = cell(3, 1);
% create datetime array with equal date intervals (1d)
% % dateIntvl =  data{:, 'date'}; % verify integrity of algorithm COMMENT
dateIntvl =  data{1, 'date'}:days(1):data{end, 'date'}; %
% n of intervals (ie days)
dateIntvlN = length(dateIntvl);
trenddata = zeros(dateIntvlN, 3);

%% Write Input Parameters to log file
writeInputLog(fID, stationname, data{:, 'date'}, ...
    polynDeg, P, HJumps, EQJump, KK, p, outl_factor);

% other variables - get count of unknowns for each type (poly, osc, jumps,
% transients)
nPolynTerms = polynDeg + 1; % 0, 1, 2, ... 
nOscParam = length(W) * 2; % cos & sin components (C, S) for every oscillation
nJumps = length(HJumps); % All Jumps - From DB and ITRF (if set to true)
nEQJumps = length(EQJump); % Only EQ Jumps -> n of transients

%% Trend Estimation
for i = 1:3
    fprintf('Evaluating "%s" ...\n', coordinateSTR{i});
    [y, result_parameterC, xEst, outlierLogical] = computeTrendIRLS(...
        t, ... % t in years where t0 = beginning of TS
        data{:, i + 2}, ... % vector with TS metric (Coordinate measurement)
        polynDeg, ...  % polynome degree
        W, ...  % periods
        HJumps, ...  % jumps: time in years since t0
        EQJump, ... % eqs jumps: time in years since t0
        T, ... %  logar. transient parameter T for earthquakes
        KK, ... % n of iterations for IRLS
        p, ... % L_p Norm for IRLS
        outl_factor); % median(error) + standard deviation * factor -> outlier
    
    % store results in master arrays for further evaluation
    trenddata(:, i) = y;
    result_parameters{i, 1} = coordinateSTR{i};
    result_parameters{i, 2} = result_parameterC; % assign parameter cell to super cell
    outlier_logicals{i} = outlierLogical; % 1: suspected outlier measurements, computed in IRLS function
    
    % print output parameters using custom function
    writeOutputLog(fID, [stationname, '-', coordinateSTR{i}], xEst, ...
        polynDeg, P, HJumps, EQJump, result_parameterC{1, 2}, result_parameterC{2, 2})
end
fclose(fID); % close log file
fprintf('Calculation finished.\nPlotting and writing results ...\n')

%% Write Trend Results to file
% use parameters in file name
resultSaveFile = fullfile(logFileFolder, [stationname, ...
    sprintf('_itrf%d_KK%d_p%.1f_outl%d', doITRFjump,KK, p, outl_factor), ...
    '.csv']); % output: file name of computed trends (csv)

% use custom fct
resultM = writeResultMatrix(dateIntvl, trenddata, doITRFjump, KK, p, outl_factor, ...
    [result_parameters{1, 2}{1, 2}, result_parameters{2, 2}{1, 2}, result_parameters{3, 2}{1, 2}], ...
    [result_parameters{1, 2}{2, 2}, result_parameters{2, 2}{2, 2}, result_parameters{3, 2}{2, 2}]);

% write matrix to csv file
%writematrix(resultM, resultSaveFile, 'Delimiter', 'comma') % R2019a
csvwrite(resultSaveFile, resultM); % R2006

%% Visualize Results
% set up title
titleString = cell(3, 1);
titleStringPattern = 'Station: %s - n of Obs.=%d - ITRF Jumps: %s - RSME = %.2fmm';

for i = 1:3 % for E N U respectively
    if doITRFjump == true
        ITRFstring = 'true';
    else
        ITRFstring = 'false';
    end
    % set up plot title
    titleString{i} = sprintf(titleStringPattern, ...
        stationname, size(data{:, 'date'}, 1), ITRFstring, result_parameters{i,2}{1,2});
end

figure
VisualizeTS_Trend_Outliers_ITRF_ENU(...
    data{:, 'date'}, [data{:, 'E'}, data{:, 'N'}, data{:, 'U'}], outlier_logicals, coordinateSTR, ...
    dateIntvl, trenddata, ...
    titleString, ...
    currStationJumps{:, 'Date'}, ...
    [currStationJumps{:, 'Earthquake'}, currStationJumps{:, 'HWSW_Change'}, currStationJumps{:, 'Unknown'}], ...
    jump_category_names, ...
    readITRFChanges(itrf_changes_textfile)...
    )

% set(gcf, 'InnerPosition', [0 0 604 513]);
set(gcf, 'InnerPosition', [0 0 1000 600]); % large figure
fprintf('Done!\n')

% close all