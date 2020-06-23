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
itrfChangesTextfile = 'src/itrf_changes.txt';
doSaveResults = false; % save pngs and result files
%%% Name of station to be analysed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTION FOR ANALYSIS

% stationName = 'MEXI'; % Mexicali, Mexico
% stationName = 'ALAR'; % Arapiraca, Brazil
% stationName = 'AREQ'; % Arequipa, Peru
% stationName = 'CONZ'; % Concepcion, Chile
% stationName = 'OAX2'; % Oaxaca, Mexico
% stationName = 'CUEC'; % Cuenca, Ecuador
% stationName = 'MZAE'; % Santa Rosa, Argentina (missing jump)
% stationName = 'NEIL'; % Ciudad Neilly, Costa Rica
% stationName = 'RWSN'; % Rawson, Argentina
% stationName = 'PBJP';

% stationName = '21701S007A03'; % KSMV %[ok]
stationName = '21702M002A07'; % MIZU %[ok]
% stationName = '21729S007A04'; % USUDA %[ok]
% stationName = '21754S001A01'; % P-Okushiri - Hokkaido %[ok, 2 eqs, doeqjumps]
% stationName = '21778S001A01'; % P-Kushiro - Hokkaido %[ok, 2 eqs, doeqjumps]
% stationName = '23104M001A01'; % Medan (North Sumatra) %[ok, 2polynDeg, 2 eqs, doeqjumps]
% stationName = '41705M003A04'; % Santiago %[ok, doeqjumps]
% stationName = '41719M004A02'; % Concepcion %[ok]

%%% Trend Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
polynDeg = [-1, -1, -1];% Polynomial Trend,Integer Degree, Range [-1..3]
% periods / oscillations in YEARS (=365.25 days) in vector form
osc = {[], [], []};     % common values: 0.5y, 1y
% Model ITRF jumps (set to "true") or ignore ITRF jumps (set to "false")
doITRFjump  = [false false false]; % E-N-U
doEQjump    = [true true true]; % E-N-U
% specify type of transient: "log","exp","nil"
transientType = {...
    'log','log'; ...    % coordinate1:E|X
    'log','log'; ...    % coordinate2:N|Y
    'log','log'};       % coordinate3:U|Z

% Parameter tau in [years] for computation of logarithmic transient for
% earthquake events (jumps):
% vector mapping different T (tau) relaxation coefficients
tauVec1 = years(days(1:10:200));
tauVec2 = years(days(250:30:730));
% tauVec1 = years(days(1:5:365)); % full range tau

% Additional Parameters for LSE/IRLSE (can be adjusted with care)
KK = 0;             % n of iterations for IRLS
p = 2.0;            % L_p Norm for IRLS
outlFactor = 100;   % median(error) + standard deviation * factor -> outlier

% other variables
coordinateName    = {'E [mm]', 'N [mm]', 'U [mm]'}; % used to label plots
% coordinateName    = {'X [mm] rel. to x_1', 'Y [mm] rel. to y_1', 'Z [mm] rel. to z_1'}; % used to label plots
jumpCategoryNames = {'Earthquake', 'SW/HW-Change', 'Unknown'}; % corresponds to jump table columns

%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logFileFolder = 'TSA_TrendComputationResults'; % output: log file directory

%% Log File
if ~exist(logFileFolder, 'dir')
    mkdir(logFileFolder)
    fprintf('Result file storage directory "%s" created.\n', logFileFolder);
end

%% Load Station Data
% look in matched files
% naming pattern for station measurement data has to be "<StationName>.mat"
fpathpattern = fullfile(inputFolder, sprintf('%s.mat', stationName)); % match pattern
if exist(fpathpattern, 'file')
    load(fpathpattern)
    fprintf('Evaluating Station "%s".\n', stationName);
else
    error('Specified station "%s" - associated .mat-file could not be found!', stationName)
end

% Reassign currStation variable
data = currStation.Data;
dataStation = currStation.Station;
data = sortrows(data, 1); % sort rows according to datetime column ASC

% Load Jump Table using custom import function
dataJump = importfileJumpCSV(jumpCSVLocation);

% Get jumps for this station
isStationJump = strcmp(dataJump{:, 'Station'}, stationName); % selection logical
fprintf('%d jumps for this station found.\n', nnz(isStationJump));

% All Jumps for current station
currStationJumps = dataJump(isStationJump, :);
% All Jumps for current station where "Use" Column Element is set to "1"
currStationJumps = currStationJumps(currStationJumps.Use == 1, :);

%% Visualise Raw Time Series and Jumps from Jump Table
% Plot TS and Jumps
figVis = figure;
% VisualizeTS_ENU(data, dataSTATION_NAME, currStationJumps{:, 2})
VisualizeTS_ENU2(data, dataStation, currStationJumps{:, 2}, ...
    'jumpTypes', currStationJumps{:, 4:5});

%% Prepare Data for LS Estimation
t  = data{:, 't'}; % [seconds]; for years, do /(365.25 * 86400);
t  = t-t(1); % adjust for negative t values
t0 = data{1, 'date'}; % beginning of ts as datetime
% convert oscillations to angular velocity
oscW = cellfun(@(x) 2*pi./x, osc, 'UniformOutput', false);

% preallocate
heavJumps  = cell(3,1); % Heaviside Jump cells for the 3 coordinates E,N,U
transients = cell(3,1); % Transient cells for the 3 coordinates E,N,U 

% load tables
% Distinguish between EQs (invokes transient) and other jumps (do not invoke transient)
eqLogical = logical(currStationJumps{:, 4}); % 1:=earthquake; 0:=no earthquake
jumps0itrf = getRelativeITRFJumps(t0, itrfChangesTextfile);

% assign tables
for i = 1:3 % E-N-U
    % Get Earthquake Transient Vector - Those Jumps will invoke a transient
    transients{i} = getRelativeJumps_eq(currStationJumps{:, 2}, t0, eqLogical);
    fprintf('[%s]: doEQjump set to "%s"\n', coordinateName{i}(1), mat2str(doEQjump(i)));
    if doEQjump(i) % all jumps
        heavJumps{i} = getRelativeJumps(currStationJumps{:, 2}, t0);
    else % jumps w/o eq
        heavJumps{i} = getRelativeJumps_eq(currStationJumps{:, 2}, t0, ~eqLogical);
    end
    % Add ITRF Realization Change Jumps
    fprintf('[%s]: doITRFjumps set to "%s"\n', coordinateName{i}(1), mat2str(doITRFjump(i)));
    if doITRFjump(i)
        % append itrf jumps to heaviside jumps vector
        heavJumps{i} = [heavJumps{i}; jumps0itrf];
    end
end

% preallocate cell arrays for output
result_parameters = cell(3, 2);
outlier_logicals  = cell(3, 1);
outlier_logicals  = cellfun(@(x) zeros(length(t), 1), outlier_logicals, 'UniformOutput', false); % workaround for no outliers

% create datetime array with equal date intervals (1d)
tInterpolV =  data{:, 'date'}; % verify integrity of algorithm COMMENT
% dateIntvl =  min(data{:, 'date'}):days(1):max(data{:, 'date'}); %
% n of intervals (ie days)
dateIntvlN = length(tInterpolV);
% trenddata = zeros(dateIntvlN, 3); % preallocate trend data array
trenddata = [];

% generate tau vector for transients containing all combinations
tsFctStr = {'log','exp'};
tauCell = cell(3,1);
tauTypes = cell(3,1);
for i = 1:3 % E-N-U
    if any(strcmp( transientType{i,1} , tsFctStr )) &&  any(strcmp( transientType{i,2} , tsFctStr ))
        [tauGrid1, ...
            TauGrid2] = meshgrid(tauVec1, tauVec2); % create two grids
        tauGrid       = cat(2, tauGrid1, TauGrid2); % cat along 2nd dimension
        tauVec = reshape(tauGrid, [], 2);           % reshape to 2 col vector with rows (tau1, tau2)
        tauCell{i} = num2cell(tauVec,2);
        tauTypes{i} = [ transientType{i,1};transientType{i,2} ];
    elseif any(strcmp( transientType{i,1} , tsFctStr )) && ~any(strcmp( transientType{i,2} , tsFctStr ))
        tauTypes{i} = [ transientType{i,1} ];
        tauVec        = tauVec1';
        tauCell{i} = num2cell(tauVec,2);
    elseif ~any(strcmp( transientType{i,1} , tsFctStr )) &&  any(strcmp( transientType{i,2} , tsFctStr ))
        tauTypes{i} = [ transientType{i,2} ];
        tauVec        = tauVec2';
        tauCell{i} = num2cell(tauVec,2);
    else % no match found, no transient to be estimated
        tauTypes{i} = '';
        tauCell{i} = {[]};
    end
%     tauCell{i} = tauVec;
end

% other variables - get count of params (poly,osc, jumps, transients) per coordinate
for i = 1:3 % E-N-U
    nPolynTerms(i) = polynDeg(i) + 1; % 0, 1, 2, ...
    nOscParam(i)   = length(oscW{i}) * 2; % cos & sin components (C, S) for every oscillation
    nJumps(i)      = length(heavJumps{i}); % All Jumps - From DB and ITRF (if set to true)
    nTransients(i) = length(transients{i}) * (size(tauCell{i}{1},2)); % Only EQ Jumps -> n of transients
end

% prepare array to store results rms,wrms,est.params FOR EVERY COMBINATION
% OF TAU (rows) and E,N,U (cols)
resultCell = cell(3,1);
for i = 1:3 % E-N-U
    resultCell{i} = zeros(size(tauCell{i},1), 2+nPolynTerms(i)+nOscParam(i)+nJumps(i)+nTransients(i) );
end

%% Write Input Parameters to log file
% for i = 1:3 % E-N-U
%     writeInputLog(fID, stationName, coordinateName{i}, data{:, 'date'}, ...
%         polynDeg(i), osc{i}, heavJumps{i}, ...
%         transients{i}, cell2mat(tauCell{i}), tauTypes{i}, ...
%         KK, p, outlFactor);
% end

%% Parameter Estimation
for j = 1:3 % E-N-U
    for i = 1:length(tauCell{j})
        % LSE to get approximate parameters x0
        fprintf('Evaluating "%s" ...\n', coordinateName{j});
        [y, result_parameterC, xEst, ~] = computeTrendIRLS(... % trend, rms/wrms, parameters, outlier logical
            t, ...                  % t in years where t0=beginning of TS
            data{:, j + 2}, ...     % vector with metric (coordinate)
            polynDeg(j), ...        % polynome degree
            oscW{j}, ...            % periods
            heavJumps{j}, ...       % jumps: time in years since t0
            transients{j}, ...      % eq transients: time in years since t0
            tauCell{j}{i}, ...      % transient parameter tau
            tauTypes{j},...         % type (function) of tau 2(log|exp)
            KK, ...                 % n of iterations for IRLS
            p, ...                  % L_p Norm for IRLS
            outlFactor);            % median(error) + standard deviation * factor -> outlier
        
        % results: [RMS,WRMS,Params]
        resultCell{j}(i,:) = [result_parameterC{1,2}, result_parameterC{1,2}, xEst'];
    end
end

% get best solution: min RMS (idx)
[~, EMinIdx]    = min( resultCell{1}(:,1) );
[~, NMinIdx]    = min( resultCell{2}(:,1) );
[~, UMinIdx]    = min( resultCell{3}(:,1) );
if length(tauCell{1})==length(tauCell{2})
    [~, ENMinIdx]   = min(sum( [resultCell{1}(:,1), resultCell{2}(:,1)],2));
    if length(tauCell{1})==length(tauCell{2}) &&  ...
            length(tauCell{1})==length(tauCell{3})
        [~, TotalMinIdx]= min(sum( [resultCell{1}(:,1), resultCell{2}(:,1), resultCell{3}(:,1)] ,2));
    end
end
BestIdx = [EMinIdx NMinIdx UMinIdx]; % can be adapted

% compute time series based on found solution for ENU
for i = 1:3 % E-N-U
    oscCS = resultCell{i}(BestIdx(i), 3+nPolynTerms(i) : 3+nPolynTerms(i)+nOscParam(i)-1);
    oscCS = [oscCS(1:2:end); oscCS(2:2:end)]; % reorder cosine/sine components
    trenddata(:, i) = TimeFunction(years(seconds(t)), ...
        resultCell{i}(BestIdx(i), 3 : 3+nPolynTerms(i)-1 ), ...
        oscCS, ... % rearranged oscillation amplitudes (cosine, sine)
        osc{i}, ... % periods
        years(seconds(heavJumps{i})), ...
        resultCell{i}(BestIdx(i), 3+nPolynTerms(i)+nOscParam(i) : 3+nPolynTerms(i)+nOscParam(i)+nJumps(i)-1), ...
        years(seconds(transients{i})), ...
        resultCell{i}(BestIdx(i), 3+nPolynTerms(i)+nOscParam(i)+nJumps(i) : end), ...
        tauCell{i}{BestIdx(i)},...
        tauTypes{i}...
        );
end

% create map plot for tau
for i = 1:3 % E-N-U
    if size(tauCell{i}{1},2) == 1 % 1 tau
        figure
        plot(days(years( cell2mat(tauCell{i}) )) , resultCell{i}(:,1) )
        xlabel('\tau_{log} [days]')
        ylabel('rms [mm]')
        title(["\tau parameter space for ", coordinateName{i}])
    elseif size(tauCell{i}{1},2) == 2 && length(tauVec1)>1 && length(tauVec2)>1  % 2 tau
        resultGrid = reshape(resultCell{i}(:,1), length(tauVec2), length(tauVec1) );
        figure
        colormap(flipud(parula)) % low error = good
        contourf(days(years(tauVec1)),days(years(tauVec2)),resultGrid);
        xlabel('\tau_{log1}SHORT [days]')
        ylabel('\tau_{log2}LONG [days]')
        title(["\tau parameter space for ", coordinateName{i}])
        c = colorbar;
        c.Label.String = 'RMS';
    end
end

% save results, write logs
for i = 1:3 % E-N-U
    result_parameters{i, 1} = coordinateName{i};
    result_parameters{i, 2} = ...
        {'rms', resultCell{i}(BestIdx(i),1); 'wrms', resultCell{i}(BestIdx(i),2)}; % assign parameter cell to super cell
end

%% Write Trend Results to files
% Open Log File and get identifier
logFile = [stationName,'.', datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS'),'.log']; % output: log file name
fID = fopen(fullfile(logFileFolder, logFile), 'wt');
for i = 1:3 % E-N-U
    writeInputLog(fID, stationName, coordinateName{i}, data{:, 'date'}, ...
        polynDeg(i), osc{i}, heavJumps{i}, ...
        transients{i}, cell2mat(tauCell{i}), tauTypes{i}, ...
        KK, p, outlFactor);
    writeOutputLog(fID, dataStation, coordinateName{i}, ...
        resultCell{i}(BestIdx(i), 3 : 3+nPolynTerms(i)-1 ), ...
        oscCS, ...
        resultCell{i}(BestIdx(i), 3+nPolynTerms(i)+nOscParam(i) : 3+nPolynTerms(i)+nOscParam(i)+nJumps(i)-1), ...
        resultCell{i}(BestIdx(i), 3+nPolynTerms(i)+nOscParam(i)+nJumps(i) : end), ...
        tauCell{i}{BestIdx(i)}, ...
        resultCell{i}(BestIdx(i),1), resultCell{i}(BestIdx(i),2)...
        )
end

fclose(fID); % close log file
fprintf('Calculation finished.\nPlotting and writing results ...\n')

% use parameters in file name
resultSaveFile = fullfile(logFileFolder, [stationName, ...
    sprintf('_itrf%d_KK%d_p%.1f_outl%d', doITRFjump,KK, p, outlFactor), ...
    '.csv']); % output: file name of computed trends (csv)

% use custom fct
resultM = writeResultMatrix(tInterpolV, trenddata, doITRFjump, KK, p, outlFactor, ...
    [result_parameters{1, 2}{1, 2}, result_parameters{2, 2}{1, 2}, result_parameters{3, 2}{1, 2}], ...
    [result_parameters{1, 2}{2, 2}, result_parameters{2, 2}{2, 2}, result_parameters{3, 2}{2, 2}]);

% write matrix to csv file
%writematrix(resultM, resultSaveFile, 'Delimiter', 'comma') % R2019a
if doSaveResults; csvwrite(resultSaveFile, resultM); end% R2006

%% Visualize Results
% set up title
titleString = cell(3, 1); % preallocate
for i = 1:3 % E-N-U
    % set up plot title
    if isempty(tauTypes{i});tsStr='none'; else; tsStr=sprintf('%s %s', transientType{i,1},transientType{i,2}); end
    titleString{i} = sprintf('Station:"%s" Transient:%s jump(itrf):%s  jump(eq):%s  RMS=%.2fmm WRMS=%.2fmm', ...
        stationName, tsStr,...
        mat2str(doITRFjump(i)), ...
        mat2str(doEQjump(i)), ...
        result_parameters{i,2}{2,2}, result_parameters{i,2}{2,2}); % rms&wrms
end

% visualize time series and results
figTSA = figure;
VisualizeTS_Trend_Outliers_ITRF_ENU(...
    data{:, 'date'}, [data{:, 3}, data{:, 4}, data{:, 5}], outlier_logicals, coordinateName, ...
    tInterpolV, trenddata, ...
    titleString, ...
    currStationJumps{:, 'Date'}, ...
    [currStationJumps{:, 'Earthquake'}, currStationJumps{:, 'HWSW_Change'}, currStationJumps{:, 'Unknown'}], ...
    jumpCategoryNames, ...
    readITRFChanges(itrfChangesTextfile)...
    )
% set(gcf, 'InnerPosition', [0 0 604 513]); % small figure
set(gcf, 'InnerPosition', [0 0 1000 600]); % large figure

% visualize residuals
figRes = figure;
VisualizeResiduals(...
    data{:, 'date'}, [...
    data{:, 3}-trenddata(:,1), ... % E residuals
    data{:, 4}-trenddata(:,2), ... % N residuals
    data{:, 5}-trenddata(:,3)], ...% U residuals
    outlier_logicals, ...
    cellfun(@(x) ['Residual \Delta',x],coordinateName,'UniformOutput',false), ...
    titleString, ...
    currStationJumps{:, 'Date'}, ...
    [currStationJumps{:, 'Earthquake'}, currStationJumps{:, 'HWSW_Change'}, currStationJumps{:, 'Unknown'}], ...
    jumpCategoryNames, ...
    readITRFChanges(itrfChangesTextfile)...
    )
% set(gcf, 'InnerPosition', [0 0 604 513]); % small figure
set(gcf, 'InnerPosition', [0 0 1000 600]); % large figure

% PRINT PLOT TREND
if doSaveResults % Save figure as image file
    plot_title = [stationName, '-trend.png'];
    plot_dir = 'stationTSA_dailyXYZfiles_xyz_plots';
    saveas(figTSA, fullfile(plot_dir, plot_title));
end

% PRINT PLOT RESIDUALS
if doSaveResults % Save figure as image file
    plot_title = [stationName, '-residuals.png'];
    plot_dir = 'stationTSA_dailyXYZfiles_xyz_plots';
    saveas(figRes, fullfile(plot_dir, plot_title));
end 

%%
fprintf('Done!\n')
% close all