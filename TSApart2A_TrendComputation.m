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
doSaveResults = false; % save pngs and result files
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

% stationname = '21701S007A03'; % KSMV %[ok]
stationname = '21702M002A07'; % MIZU %[ok]
% stationname = '21729S007A04'; % USUDA %[ok]
% stationname = '21754S001A01'; % P-Okushiri - Hokkaido %[ok, 2 eqs, doeqjumps]
% stationname = '21778S001A01'; % P-Kushiro - Hokkaido %[ok, 2 eqs, doeqjumps]
% stationname = '23104M001A01'; % Medan (North Sumatra) %[ok, 2polynDeg, 2 eqs, doeqjumps]
% stationname = '41705M003A04'; % Santiago %[ok, doeqjumps]
% stationname = '41719M004A02'; % Concepcion %[ok]

%%% Trend Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial Trend: Degree
polynDeg = -1; % integer degree number
% polynDeg = 2;
% polynDeg = 3;

% periods / oscillations in YEARS (=365.25 days) in vector form
P = [];
% P(1) = 1;
% P(2) = 1/2;
% P(3) = 10;

% convert oscillations to angular velocity
W = 2*pi./P;

% Parameter T in [years] for computation of logarithmic transient for
% earthquake events (jumps)
T = years(days(10));
% vector mapping different T (tau) relaxation coefficients
tauVec1 = years(days(1:10:120));
tauVec2 = years(days(201:30:730));
% tauVec2=[]; % only 1 transient

% Model ITRF jumps (set to "true") or ignore ITRF jumps (set to "false")
doITRFjump = false; % E - N - U
doEQjump = true;

% Additional Parameters for LSE/IRLSE (can be adjusted with care)
KK = 0; % n of iterations for IRLS
p = 2.0; % L_p Norm for IRLS
outl_factor = 100; % median(error) + standard deviation * factor -> outlier

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
coordinateSTR = {'E [mm]', 'N [mm]', 'U [mm]'}; % used to label plots
% coordinateSTR = {'X [mm] rel. to x_1', 'Y [mm] rel. to y_1', 'Z [mm] rel. to z_1'}; % used to label plots

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
data = sortrows(data, 1); % sort rows according to datetime column ASC

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
t=t-t(1); % adjust for negative t values
t0 = data{1, 'date'}; % beginning of ts as datetime

%% Jump table(s)
% Distinguish between EQs (invokes log. transient) and other jumps (unknown cause or HW
% change -> do not invoke transient)
EQLogical = logical(currStationJumps{:, 4}); % 1:=earthquake; ~1:=no earthquake
% Get Earthquake Jump Vector - Those Jumps will invoke a log. transient
EQJump = getRelativeJumps_eq(currStationJumps{:, 2}, t0, EQLogical);

if doEQjump
    HJumps = getRelativeJumps(currStationJumps{:, 2}, t0);
else
    HJumps = getRelativeJumps_eq(currStationJumps{:, 2}, t0, ~EQLogical);
end
%% Add ITRF Realization Change Jumps
fprintf('doITRFjumps set to "%s"\n',  doITRFjump);
if doITRFjump
    jumps0itrf = getRelativeITRFJumps(t0, itrf_changes_textfile);
    % append itrf jumps to heaviside jumps vector
    HJumps = [HJumps; jumps0itrf];
end

%% Prepare Trend Estimation
% preallocate cell arrays
result_parameters = cell(3, 2);
outlier_logicals = cell(3, 1);
outlier_logicals = cellfun(@(x) zeros(length(t), 1), outlier_logicals, 'UniformOutput', false); % workaround for no outliers

% create datetime array with equal date intervals (1d)
tInterpolV =  data{:, 'date'}; % verify integrity of algorithm COMMENT
% dateIntvl =  min(data{:, 'date'}):days(1):max(data{:, 'date'}); %
% n of intervals (ie days)
dateIntvlN = length(tInterpolV);
% trenddata = zeros(dateIntvlN, 3); % preallocate trend data array
trenddata = [];

%% Write Input Parameters to log file
writeInputLog(fID, stationname, data{:, 'date'}, ...
    polynDeg, P, HJumps, EQJump, KK, p, outl_factor);

% other variables - get count of unknowns for each type (poly, osc, jumps,
% transients)
nPolynTerms = polynDeg + 1; % 0, 1, 2, ...
nOscParam = length(W) * 2; % cos & sin components (C, S) for every oscillation
nJumps = length(HJumps); % All Jumps - From DB and ITRF (if set to true)
nEQJumps = length(EQJump) * (size(tauVec1, 1)+size(tauVec2, 1)); % Only EQ Jumps -> n of transients
%% Preparation for Trend Estimation
% set up tau vector
if ~isempty(tauVec2) % only if two parameters per EQ event
    [tauGrid1, TauGrid2] = meshgrid(tauVec1, tauVec2); % create two grids
    tauGrid = cat(2, tauGrid1, TauGrid2); % cat along 2nd dimension
    tauVec = reshape(tauGrid, [], 2); % reshape to 2 col vector with rows (tau1, tau2)
else % if only 1 parameter tau per EQ event
    tauVec = tauVec1';
end
% prepare array to store results rms,wrms,est.params FOR EVERY COMBINATION
% OF TAU and E,N,U
resultArray = zeros( size(tauVec,2) , 2+nPolynTerms+nOscParam+nJumps+nEQJumps , 3 );
%% Trend Estimation
% Loop tau values
for i = 1:length(tauVec)
    for j = 1:3 % E,N,U
        % LSE to get approximate parameters x0
        fprintf('Evaluating "%s" ...\n', coordinateSTR{j});
        [y, result_parameterC, xEst, ~] = computeTrendIRLS(... % trend, rms/wrms, parameters, outlier logical
            t, ... % t in years where t0 = beginning of TS
            data{:, j + 2}, ... % vector with TS metric (Coordinate measurement)
            polynDeg, ...  % polynome degree
            W, ...  % periods
            HJumps, ...  % jumps: time in years since t0
            EQJump, ... % eqs jumps: time in years since t0
            tauVec(i, :), ... %  logar. transient parameter T for earthquakes
            KK, ... % n of iterations for IRLS
            p, ... % L_p Norm for IRLS
            outl_factor); % median(error) + standard deviation * factor -> outlier
        
        resultArray(i, :, j) = [result_parameterC{1,2}, result_parameterC{1,2}, xEst'];        
    end
end
% create map plot
for i = 1:3
    if isempty(tauVec2)
        figure
        plot(days(years(tauVec)), resultArray(:,1,i))
        xlabel('\tau_{log} [days]')
        ylabel('rms [mm]')
        title(coordinateSTR{i})
    elseif ~isempty(tauVec2)
        resultGrid = reshape(resultArray(:, 1, i), length(tauVec2), length(tauVec1));
        figure
        colormap(flipud(parula)) % low error = good
        contourf(days(years(tauVec1)),days(years(tauVec2)),resultGrid);
        xlabel('\tau_{log1}SHORT [days]')
        ylabel('\tau_{log2}LONG [days]')
        colorbar
    end
end
% get best solution vTv (idx)
[~, EMinIdx]    = min( resultArray(:,1,1) );
[~, NMinIdx]    = min( resultArray(:,1,2) );
[~, UMinIdx]    = min( resultArray(:,1,3) );
[~, ENMinIdx]   = min(sum( resultArray(:, 1, 1:2) ,3));
[~, TotalMinIdx]= min(sum( resultArray(:, 1, 1:3) ,3));

% compute time series based on found solution for ENU
for i = 1:3
    trenddata(:, i) = TimeFunction(years(seconds(t)), ...
        resultArray(ENMinIdx,3 : 3+nPolynTerms-1, i), ...
        [], [], ...
        years(seconds(HJumps)), ...
        resultArray(ENMinIdx,3+nPolynTerms+nOscParam : 3+nPolynTerms+nOscParam+nJumps-1, i), ...
        years(seconds(EQJump)), ...
        resultArray(ENMinIdx,3+nPolynTerms+nOscParam+nJumps : end, i), ...
        tauVec(ENMinIdx, :));
end

% save results, write logs
for i = 1:3
    result_parameters{i, 1} = coordinateSTR{i};
    result_parameters{i, 2} = ...
        {'rms', resultArray(ENMinIdx,1,i); 'wrms', resultArray(ENMinIdx,2,i)}; % assign parameter cell to super cell
end
% writeOutputLog(fID, [stationname, '-', coordinateSTR{j}], xEst, ...
%     polynDeg, P, HJumps, EQJump, results{1, 2}, results{2, 2})
% outlier_logicals{j} = outlierLogical; % 1: suspected outlier measurements, computed in IRLS function

% continue evaluation of found solution

fclose(fID); % close log file
fprintf('Calculation finished.\nPlotting and writing results ...\n')

%% Write Trend Results to file
% use parameters in file name
resultSaveFile = fullfile(logFileFolder, [stationname, ...
    sprintf('_itrf%d_KK%d_p%.1f_outl%d', doITRFjump,KK, p, outl_factor), ...
    '.csv']); % output: file name of computed trends (csv)

% use custom fct
resultM = writeResultMatrix(tInterpolV, trenddata, doITRFjump, KK, p, outl_factor, ...
    [result_parameters{1, 2}{1, 2}, result_parameters{2, 2}{1, 2}, result_parameters{3, 2}{1, 2}], ...
    [result_parameters{1, 2}{2, 2}, result_parameters{2, 2}{2, 2}, result_parameters{3, 2}{2, 2}]);

% write matrix to csv file
%writematrix(resultM, resultSaveFile, 'Delimiter', 'comma') % R2019a
if doSaveResults; csvwrite(resultSaveFile, resultM); end% R2006

%% Visualize Results
% set up title
titleString = cell(3, 1);
titleStringPattern = 'Station: %s - n of Obs.=%d - ITRF Jumps: %s - RMS = %.2fmm';

for i = 1:3 % for E N U respectively
    if doITRFjump == true; ITRFstring = 'true';
    else; ITRFstring = 'false'; end
    % set up plot title
    titleString{i} = sprintf(titleStringPattern, ...
        stationname, size(data{:, 'date'}, 1), ITRFstring, result_parameters{i,2}{1,2});
end

figTSA = figure;
VisualizeTS_Trend_Outliers_ITRF_ENU(...
    data{:, 'date'}, [data{:, 3}, data{:, 4}, data{:, 5}], outlier_logicals, coordinateSTR, ...
    tInterpolV, trenddata, ...
    titleString, ...
    currStationJumps{:, 'Date'}, ...
    [currStationJumps{:, 'Earthquake'}, currStationJumps{:, 'HWSW_Change'}, currStationJumps{:, 'Unknown'}], ...
    jump_category_names, ...
    readITRFChanges(itrf_changes_textfile)...
    )

% set(gcf, 'InnerPosition', [0 0 604 513]);
set(gcf, 'InnerPosition', [0 0 1000 600]); % large figure

%% PRINT PLOT
plot_title = [stationname, '-trend.png'];
plot_dir = 'stationTSA_dailyXYZfiles_xyz_plots';
if doSaveResults; saveas(figTSA, fullfile(plot_dir, plot_title)); end% Save figure as image file

%%
fprintf('Done!\n')

% close all