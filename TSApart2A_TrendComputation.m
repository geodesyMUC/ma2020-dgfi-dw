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
%   Disclaimer: Certain parameters for the (Iteratively Reweighted) Least
%   Squares can be adjusted in the corresponding function file
%   ("computeTrendIRLS.mat")

% David Wallinger, DGFI, 5.8.2019

clear variables;
close all;
addpath('myfunctions')

%% SETTINGS (adapt if necessary) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INPUT data settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputFolder = 'station_data'; % Where Station Data (TSA_ReadAndTransform) is stored as ".mat"
jumpCSVLocation = 'jumps_version3.csv'; % Location of Jump Table/Jump Database

%%% Name of station to be analysed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTION FOR ANALYSIS

% stationname = 'MEXI'; % Mexicali, Mexico
% stationname = 'ALAR'; % Arapiraca, Brazil
% stationname = 'AREQ'; % Arequipa, Peru
% stationname = 'CONZ'; % Concepcion, Chile
% stationname = 'OAX2'; % Oaxaca, Mexico
% stationname = 'CUEC'; % Cuenca, Ecuador
stationname = 'MZAE'; % Santa Rosa, Argentina (missing jump)
% stationname = 'NEIL'; % Ciudad Neilly, Costa Rica
% stationname = 'RWSN'; % Rawson, Argentina
% stationname = 'PBJP'; %

%%% Trend Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial Trend: Degree
polynDeg = 1; % integer degree number
% polynDeg = 2;
% polynDeg = 3;

% periods / oscillations in YEARS (=365.25 days) in vector form
P = [];
P(1) = 1;
P(2) = 1/2;
% P(3) = 10;

% Parameter T in [years] for computation of logarithmic transient for
% earthquake events (jumps)
T = 1;

% Model ITRF jumps (set to "true") or ignore ITRF jumps (set to "false")
doITRFjump = true; % E - N - U

% Additional Parameters for LSE/IRLSE (can be adjusted with care)
KK = 10; % n of iterations for IRLS
p = 1.5; % L_p Norm for IRLS
outl_factor = 5; % median(error) + standard deviation * factor -> outlier

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
% reason: Numerical Stability in LSE matrix inversion
t = data{:, 't'}; % [seconds]; for years, do /(365.25 * 86400);
t0 = data{1, 'date'}; % beginning of ts as datetime

% convert oscillations to angular velocity
W = 2 * pi ./ P;

%% Jump table
HJumps = getRelativeJumps(currStationJumps{:, 2}, t0);

% Distinguish between EQs (invokes log. transient) and other jumps (unknown cause or HW
% change -> do not invoke transient)
EQLogical = logical(currStationJumps{:, 4}); % 1:=earthquake; ~1:=no earthquake
% Get Earthquake Jump Vector - Those Jumps will invoke a log. transient
EQJump = getRelativeJumps_eq(currStationJumps{:, 2}, t0, EQLogical); 

%% Add ITRF Realization Change Jumps
% format yyyy-MM-dd, 
ditrf(1) = datetime('2011-04-17', 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
ditrf(2) = datetime('2017-01-29', 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
% can be extended! keep format "yyyy-MM-dd", time will be 0:00am)

if doITRFjump   
    fprintf('Considering ITRF jumps, doITRFjump is set to "true".\n')
    jumps0itrf = etime(datevec(ditrf), ...
        datevec(repmat(t0, size(ditrf, 1), 1)));
    % remove negative values -> itrf jumps BEFORE Time Series starts
    jumps0itrf = jumps0itrf(jumps0itrf >= 0);
    % convert jump times to years (365.25!!!! days) and sort asc
    jumps0itrf = sort(jumps0itrf);% ./(365.25 * 86400));
    
    % append to heaviside jumps vector
    HJumps = [HJumps; jumps0itrf]; % HeavisideJump
else
    fprintf('Not considering ITRF jumps, doITRFjump is set to "false".\n')
end

%% 



%% Prepare Trend Estimation
% preallocate arrays ---
result_parameters = cell(3, 2);
OUTLIERLOGICAL = cell(3, 1);
% create array with equal date intervals (1d)
% % dateIntvl =  data{:, 'date'}; % verify integrity of algorithm COMMENT
dateIntvl =  data{1, 'date'}:days(1):data{end, 'date'}; %
% n of intervals (ie days)
dateIntvlN = length(dateIntvl);
trenddata = zeros(dateIntvlN, 3);
%---

% Write Input Parameters to log file
fprintf(fID, 'TSA: Trend for Station "%s" -----------\n', stationname);
fprintf(fID, 'Time Series Start t0  = %s\n\n', datestr(t0, 'yyyy-mm-dd HH:MM'));

fprintf(fID, '\n### INPUT PARAMETERS ###\n');

fprintf(fID, 'Polynomial Trend Model (SLTM):\nDegree = %d\n\n', polynDeg) ;
fprintf(fID, 'Oscillations:\n%s\n\n', ...
    sprintf('w(%d) = %.1f y | %.2f d\n', [(1:length(P)); P; P.*365.25]));
% All jumps
fprintf(fID, 'Jump table (Heaviside):\n');
for i = 1:length(HJumps)
    % Print Jump Datetime Information
    fprintf(fID, 'J(%d) = t0 + %.2f y | %s\n', i, HJumps(i), ...
        datestr(t0 + seconds(HJumps(i)), 'yyyy-mm-dd HH:MM'));
end

% EQ
fprintf(fID, '\nEarthquake Jump table (Heaviside + Logarithmic Transient):\n');
for i = 1:length(EQJump)
    % Print Jump Datetime Information
    fprintf(fID, 'J(%d) = t0 + %.2f y | %s\n', i, EQJump(i), ...
        datestr(t0 + seconds(EQJump(i)), 'yyyy-mm-dd HH:MM'));
end
% LSE/IRLSE
fprintf(fID, ['\nIRLSE/LSE Parameters:\nKK = %d (n of Iterations in IRLSE)\np = %.1f (L_p Norm used for IRLS)\n', ...
    'outl_factor = %d (median of error + standard deviation * factor < outlier)\n'], KK, p, outl_factor);


fprintf(fID, '\n### LEAST SQUARE ESTIMATION ###\n');
fprintf(fID, 'Using Moore-Penrose pseudoinverse for Inversion of normal equation matrix when needed\n');
% other variables - get count of unknowns for each type (poly, osc, jumps,
% transients)
nPolynTerms = polynDeg + 1; % 0, 1, 2, ... 
nOscParam = length(W) * 2; % cos & sin components (C, S) for every oscillation
nJumps = length(HJumps); % All Jumps - From DB and ITRF (if set to true)
nEQJumps = length(EQJump); % Only EQ Jumps -> n of transients

%% Trend Estimation INCLUDING EQ Transients
fprintf(fID, '\n### OUTPUT RESULTS ###\n');
fprintf(fID, '-- Logarithmic Transients for Earthquakes estimated --\n\n');

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
    OUTLIERLOGICAL{i} = outlierLogical; % 1: suspected outlier measurements, computed in IRLS function
    
    % print stuff
    strP = sprintf('p(%d) = %.5f | ', [1:nPolynTerms; xEst(1:nPolynTerms)']); % Polynomial Coefficients
    if ~isempty(P)
        strW = sprintf('w(%d): A = % .2fmm, C = % .2fmm,  S = % .2fmm\n', ...
            [1:length(P); ...
            sqrt(xEst(nPolynTerms + 1:2:nPolynTerms + nOscParam)'.^2 + ...
            xEst(nPolynTerms + 2:2:nPolynTerms + nOscParam)'.^2); % Amplitude A
            xEst(nPolynTerms + 1:2:nPolynTerms + nOscParam)'; ... % Oscillations C
            xEst(nPolynTerms + 2:2:nPolynTerms + nOscParam)']); % Oscillations S
    else
        strW = 'None estimated';
    end
    if ~isempty(HJumps)
        strH = sprintf('J(%d) = % .2fmm\n', [1:nJumps; ...
            xEst(nPolynTerms + nOscParam + 1:nPolynTerms + nOscParam + nJumps)']); % Jumps
    else
        strH = 'None estimated';
    end
    if ~isempty(EQJump)
        strEQ = sprintf('A(%d) = % .2fmm\n', [1:nEQJumps; ...
            xEst(...
            nPolynTerms + nOscParam + nJumps + 1:...
            nPolynTerms + nOscParam + nJumps + nEQJumps)']); % EQ Jumps
    else
        strEQ = 'None estimated';
    end
    % Print all Parameters to log file
    fprintf(fID, ['Estimated Parameters for %s:\nPolynomial Coefficients:\n%s\n', ...
        'Oscillation Amplitudes:\n%s\nHeaviside Jumps:\n%s\nLogarithmic Transients:\n%s\nRoot Mean Square Error RSME = %.3fmm\n\n'], ...
        coordinateSTR{i}, ...
        strP, strW, strH, strEQ, result_parameterC{1, 2});
end

fprintf('Calculation finished.\nPlotting and writing results ...\n')

%% Write Trend Results to file
% use parameters in file name
resultSaveFile = fullfile(logFileFolder, [stationname, ...
    sprintf('_itrf%d_KK%d_p%.1f_outl%d', doITRFjump,KK, p, outl_factor), ...
    '.csv']); % output: file name of computed trends (csv)

% Set up matrix with results and LSE parameters
resultM = [posixtime(dateIntvl'),  trenddata];
resultM = [[200, result_parameters{1, 2}{2, 2}, result_parameters{2, 2}{2, 2}, result_parameters{3, 2}{2, 2}]; ...
    resultM]; % 200: WRMS (third line)
resultM = [[100, result_parameters{1, 2}{1, 2}, result_parameters{2, 2}{1, 2}, result_parameters{3, 2}{1, 2}]; ...
    resultM]; % 100: RMS (second line)
resultM = [[doITRFjump, KK, p, outl_factor]; resultM]; % LSE parameters (first line)

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

% Plot Time Series & Trends
figure
VisualizeTS_Trend_Outliers_ITRF_ENU(data{:, 'date'}, ...
    [data{:, 'E'}, data{:, 'N'}, data{:, 'U'}], ...
    dateIntvl, trenddata, ...
    titleString, ...
    OUTLIERLOGICAL, ...
    currStationJumps, ...
    ditrf);
% set(gcf, 'InnerPosition', [0 0 604 513]);
set(gcf, 'InnerPosition', [0 0 1000 600]); % large figure

% close log file
fclose(fID);
fprintf('Done!\n')

% close all