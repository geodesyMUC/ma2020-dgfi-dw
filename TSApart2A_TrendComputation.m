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
doStaticFile = true;
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
% stationName = '21702M002A07'; % MIZU %[ok]
% stationName = '21729S007A04'; % USUDA %[ok]
% stationName = '21754S001A01'; % P-Okushiri - Hokkaido %[ok, 2 eqs, doeqjumps]
% stationName = '21778S001A01'; % P-Kushiro - Hokkaido %[ok, 2 eqs, doeqjumps]
stationName = '23104M001A01'; % Medan (North Sumatra) %[ok, 2polynDeg, 2 eqs, doeqjumps]
% stationName = '41705M003A04'; % Santiago %[ok, doeqjumps]
% stationName = '41719M004A02'; % Concepcion %[ok]

%%% Trend Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
polynDeg = [1, 1, 1];% Polynomial Trend,Integer Degree, Range [-1..3]
% periods / oscillations in YEARS (=365.25 days) in vector form
osc = {[], [], []};     % common values: 0.5y, 1y
% Model ITRF jumps (set to "true") or ignore ITRF jumps (set to "false")
doITRFjump  = [false false false]; % E-N-U
doEQjump    = [true true true]; % E-N-U
doStopTs = true; % global

% specify type of transient: "log","exp","nil"
transientType = {...
    'log','log'; ...    % coordinate1:E|X
    'log','log'; ...    % coordinate2:N|Y
    'log','log'};       % coordinate3:U|Z

% Parameter tau in [years] for computation of logarithmic transient for
% earthquake events (jumps):
% vector mapping different T (tau) relaxation coefficients [years]
tauVec1 = years(days(1:20:200));
tauVec2 = years(days(250:40:730));
% tauVec1 = years(days(1:5:365)); % full range tau

% optimization constraints for transient 1 and transient 2
lowLimit = [min(tauVec1), min(tauVec2)];
uppLimit = [max(tauVec1), max(tauVec2)];

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
% Convert Oscillation Periods to Angular Velocity
oscW = cellfun(@(x) 2*pi./x, osc, 'UniformOutput', false);

% Time Series Timestamps
t  = data{:, 't'};       % timestamps [seconds];
t = years( seconds(t) ); % convert timestamps to [years]
t  = t-t(1);             % adjust for negative t values
t0 = data{1, 'date'};    % timestamp [datetime] of first measurement

% Set up variables for relative date of jumps and transients
% Distinguish between EQs (invokes transient) and other jumps (do not invoke transient)
isEq = logical(currStationJumps{:, 4}); % 1:=earthquake; 0:=no earthquake
jumps0itrf = getRelativeITRFJumps(t0, itrfChangesTextfile);

% Preallocate Storages
heavJumps  = cell(3,1); % Heaviside Jumps for the 3 coordinates E,N,U
transients = cell(3,1); % Transients for the 3 coordinates E,N,U
for i = 1:3 % E-N-U
    % Get Transient Datetime from jumps where type = eq
    ts_t = getRelativeJumps_eq(currStationJumps{:, 2}, t0, isEq);
    fprintf('[%s]: doEQjump set to "%s"\n', coordinateName{i}(1), mat2str(doEQjump(i)));
    if doEQjump(i) 
        % all jump types
        jumps = getRelativeJumps(currStationJumps{:, 'Date'}, t0);
    else
        % only jump types which are NOT eq
        jumps = getRelativeJumps_eq(currStationJumps{:, 2}, t0, ~isEq);
    end
    
    % Add ITRF Realization Change Jumps
    fprintf('[%s]: doITRFjumps set to "%s"\n', coordinateName{i}(1), mat2str(doITRFjump(i)));
    if doITRFjump(i)
        % append itrf jumps to heaviside jumps vector
        jumps = [jumps; jumps0itrf];
    end
    
    % Convert to [years] & assign to cell array 
    transients{i} = years( ts_t );
    heavJumps{i}  = years( jumps );
    % Lookup Table for relative transient timestamps, types and limits (->optimization)
    tsLUT{i} = getTransientReferences(transients{i}, transientType(i,:), ...
        lowLimit, uppLimit, doStopTs);
end

% create datetime array with equal date intervals (1d)
tInterpol =  data{:, 'date'}; % take original timestamps, no interpolation
% tInterpol =  min(data{:, 'date'}):days(1):max(data{:, 'date'}); % 1d interpolation
dateIntvlN = length(tInterpol); % n of timestamps


% preallocate output storages
resStor = cell(1,3); % E-N-U
nMethod = 3; % gs,dhs,ip
fields = {'optp', 'xmin', 'model', 'trend', 'error', 'outl'};
vals = {cell(1,nMethod), cell(1,nMethod), cell(1,nMethod), cell(1,nMethod), cell(1,nMethod), cell(1,nMethod)};
for i = 1:3 % E-N-U
    resStor{i} = struct(...
        fields{1},vals{1}, ...
        fields{2},vals{2}, ...
        fields{3},vals{3}, ...
        fields{4},vals{4}, ...
        fields{5},vals{5}, ...
        fields{6},vals{6} ...
        );
end

result_parameters = cell(3, 2);
outlier_logicals  = cell(3, 1);
outlier_logicals  = cellfun(@(x) zeros(length(t), 1), outlier_logicals, 'UniformOutput', false); % workaround for no outliers

% generate tau vector for transients containing all combinations of values
tsFctStr = {'log','exp'}; % used to verify user input string
tauCell = cell(3,1);
tauTypes = cell(3,1);
for i = 1:3 % E-N-U
    if any(strcmp( transientType{i,1} , tsFctStr )) &&  any(strcmp( transientType{i,2} , tsFctStr ))
        % Two Transients
        [tauGrid1, tauGrid2] ...
                    = meshgrid(tauVec1, tauVec2); % create two grids
        tauGrid     = cat(2, tauGrid1, tauGrid2); % cat along 2nd dimension
        tauVec      = reshape(tauGrid, [], 2);    % reshape to 2 col vector with rows (tau1, tau2)
        tauCell{i}  = num2cell(tauVec,2);
        tauTypes{i} = {transientType{i,1}, transientType{i,2}};
        nTau(i)     = 2;
    elseif any(strcmp( transientType{i,1} , tsFctStr )) && ~any(strcmp( transientType{i,2} , tsFctStr ))
        % One Transient (tau_short)
        tauTypes{i} = transientType{i,1};
        tauVec      = tauVec1';
        tauCell{i}  = num2cell(tauVec,2);
        nTau(i)     = 1;
    elseif ~any(strcmp( transientType{i,1} , tsFctStr )) &&  any(strcmp( transientType{i,2} , tsFctStr ))
        % One Transient (tau_long)
        tauTypes{i} = transientType{i,2};
        tauVec      = tauVec2';
        tauCell{i}  = num2cell(tauVec,2);
        nTau(i)     = 1;
    else
        % No Transients, no match found
        tauTypes{i} = '';
        tauCell{i} = {[]};
        nTau(i) = 0;
    end
%     tauCell{i} = tauVec;
end

% other variables - get count of params (poly,osc, jumps, transients) per coordinate
for i = 1:3 % E-N-U
    nPolynTerms(i) = polynDeg(i) + 1; % 0, 1, 2, ...
    nOscParam(i)   = length(oscW{i}) * 2; % cos & sin components (C, S) for every oscillation
    nJumps(i)      = length(heavJumps{i}); % All Jumps - From DB and ITRF (if set to true)
    nEq(i)         = length(transients{i}); % n of eq events
    nTransients(i) = nEq(i) * nTau(i); % Only EQ Jumps -> n of transients
end

% prepare array to store results rms,wrms,est.params FOR EVERY COMBINATION
% OF TAU (rows) and E,N,U (cols)
fxResAll = cell(3,1);
for i = 1:3 % E-N-U
    fxResAll{i} = zeros( size(tauCell{i},1) , 1 );
end

%% Parameter Estimation (Grid Search)
params.t        = t;                    % timestamps [years], must be sorted
params.kk       = KK;                   % n of iterations for IRLS
params.p        = p;                    % L_p Norm for IRLS
params.outl     = outlFactor;           % median(error) + standard deviation * factor -> outlier
for i = 1:3 % E-N-U
    params.b        = data{:,i+2};      % vector with metric (coordinate)
    params.poly     = polynDeg(i);      % integer polynome degree
    params.w        = oscW{i};          % periods [angular velocity]
    params.jt       = heavJumps{i};     % jump timestamps [years] relative to t0
    
    % The following parameters need to be repeated for each transient or eq
    % to satisfy LS Input req.
    params.tst      = ...
        repelem( transients{i}' , nTau(i) );     % eq transients: time in years since t0
    params.tstype   = ...
        repmat( tauTypes{i} , [1,nEq(i)] );      % type (function) of tau (log|exp)
    
    fxRes = zeros( size(tauCell{i},1) , 1 ); % preallocate
    
    % tau value combination loop (grid search)
    for j = 1:length(tauCell{i})
        params.tau = ...
            repmat( tauCell{i}{j} , [1,nEq(i)] ); % transient parameter tau
        
        fprintf('LS for "%s" (%d of %d), tau=[%s\b] ...\n', ...
            coordinateName{i}, j, length(tauCell{i}), ...
            sprintf('%.4f ', tauCell{i}{j} ));
        
        [~, result_parameterC, ~, ~] = computeTrendIRLS( params ); % LS
        fxRes(j) = result_parameterC{1,2};       % 1=RMS, 2=WRMS
        fxResAll{i}(j) = result_parameterC{1,2}; % 1=RMS, 2=WRMS, for map plot
    end
    [~, minIdx] = min( fxRes ); % minimization min( f(x) )
    minIdx = minIdx(1);         % prevent duplicates
    
    params.tau = repmat( tauCell{i}{minIdx} , [1,nEq(i)] );       % optimum tau
    [~, result_parameterC, xEst, ~] = computeTrendIRLS( params ); % LS: get parameters
    
    optP = getParamStruct([result_parameterC{1,2}, result_parameterC{2,2}, xEst'],...
        2, nPolynTerms(i), nOscParam(i), nJumps(i), nTransients(i));
    
    trenddata = TimeFunction(... % calculate points
        t, ...
        optP.polyn, ...
        reorderSinCos( optP.oscil ), ... % rearranged oscillation amplitudes (cosine, sine)
        osc{i}, ...                      % periods
        heavJumps{i}, ...
        optP.jumps, ...
        repelem( transients{i}' , nTau(i) ), ...
        optP.tsamp, ...
        repmat( tauCell{i}{minIdx} , [1,nEq(i)] ), ...
        repmat( tauTypes{i} , [1,nEq(i)] )...
        );
    e = {'rms', optP.error(1); 'wrms', optP.error(2)};
    
    resStor{i}(1).trend = trenddata;
    resStor{i}(1).error = e;
    resStor{i}(1).optp  = optP;
    resStor{i}(1).xmin  = repmat( tauCell{i}{minIdx} , [1,nEq(i)] );
    resStor{i}(1).model = params;
    
end

%% Mapping&Parameter Optimization (DHS/IP)
restartScale = 0.01; %
tol = 0.001;
for i = 1:3 % E-N-U
    tauArray = cell2mat(tauCell{i});
    lLim{i}  = tsLUT{i}{:,'lBound'};
    uLim{i}  = tsLUT{i}{:,'uBound'};
    x0{i}    = tsLUT{i}{:,'lBound'}+1e-3;
    steps{i} = tsLUT{i}{:,'uBound'}-tsLUT{i}{:,'lBound'}; % if max(tauArray)==inf?
end

% create map plot for tau & optimization
for i = 1:3 % E-N-U
    if size(tauCell{i}{1},2) == 1 % 1 tau
        figure
        plot(days(years( cell2mat(tauCell{i}) )) , resultCell{i}(:,1) )
        xlabel('\tau_{log} [days]')
        ylabel('rms [mm]')
        title(["\tau parameter space for ", coordinateName{i}])
    elseif size(tauCell{i}{1},2) == 2 && length(tauVec1)>1 && length(tauVec2)>1  % 2 tau
        resultGrid = reshape(fxResAll{i}(:,1), length(tauVec2), length(tauVec1) );
        figure
        colormap(flipud(parula)) % low error = good
        % contourf(days(years(tauVec1)),days(years(tauVec2)),resultGrid);
        contourf(tauVec1,tauVec2,resultGrid);
        xlabel('\tau_{log1}SHORT [y]')
        ylabel('\tau_{log2}LONG [y]')
        title(["\tau parameter space for ", coordinateName{i}])
        c = colorbar;
        c.Label.String = 'RMS[mm]';
        xlim([min(tauVec1) , max(tauVec1)])
        ylim([min(tauVec2) , max(tauVec2)])
        hold off
    end
    
    % function to calculate points
    timeFn = @(a, b, c, d, e) TimeFunction(...
        t, ...
        a, ...                  % polynome term coefficients
        reorderSinCos( b ), ... % rearranged oscillation amplitudes (cosine, sine)
        osc{i}, ...             
        heavJumps{i}, ...
        c, ...                  % heaviside jumps
        tsLUT{i}{:,'time'}, ...
        d, ...                  % amplitudes
        e,...                   % tau
        tsLUT{i}{:,'type'}...
        );
    
    % optimization default parameters
    params.b        = data{:,i+2};  % vector with metric (coordinate)
    params.poly     = polynDeg(i);  % polynome degree
    params.w        = oscW{i};      % periods
    params.jt       = heavJumps{i}; % jumps: time in years relative to t0
    params.tst      = tsLUT{i}{:,'time'}';% eq transients: time in years since t0
    params.tstype   = tsLUT{i}{:,'type'}';  % type (function) of tau (log|exp)
    
    % function for optimization
    optFun = @(x) getTrendError(...
        params.t, ...     % t in years where t0=beginning of TS
        params.b, ...     % vector with metric (coordinate)
        params.poly, ...  % polynome degree
        params.w, ...     % periods
        params.jt, ...    % jumps: time in years since t0
        params.tst, ...   % eq transients: time in years since t0
        x, ...            % transient parameter tau
        params.tstype,... % type (function) of tau 2(log|exp)
        params.kk, ...    % n of iterations for IRLS
        params.p, ...     % L_p Norm for IRLS
        params.outl...    % median(error) + standard deviation * factor -> outlier );
        );
    
    for j = 2:3 % dhs-ip
        if j == 2
            % DHS (custom method)
            [xMin, fxMin(i),nFnCalls(i),nRestarts(i),~] = ...
                dhscopt(optFun, x0{i}', steps{i}', lLim{i}', uLim{i}', tol, restartScale, false);
        elseif j == 3
            % MATLAB method
            if ~isempty(x0{i})
                options = optimoptions(@fmincon,...
                    'Display','iter','Algorithm','interior-point');
                [xMin,fxMin_{i}] = fmincon(optFun,x0{i},...
                    [],[],[],[],lLim{i},uLim{i},[],options);
            else
                xMin = [];
            end
        end
        
        params.tau = xMin; % set transient parameters tau
        [~, result_parameterC, xEst, ~] = computeTrendIRLS( params ); % get parameters
        optP = getParamStruct([result_parameterC{1,2}, result_parameterC{2,2}, xEst'], ...
            2, nPolynTerms(i), nOscParam(i), nJumps(i), length(xMin));
        optY = timeFn(...
            optP.polyn, ...
            optP.oscil, ...
            optP.jumps, ...
            optP.tsamp, ...
            xMin);
        e = {'rms', optP.error(1); 'wrms', optP.error(2)};
        % store in result struct
        resStor{i}(j).optp = optP;
        resStor{i}(j).error = e;
        resStor{i}(j).trend = optY;
        resStor{i}(j).xmin = xMin;
        resStor{i}(j).model = params;
    end
end

%% Evaluation
if ~doStaticFile
    currTime = [datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS'), '.'];
else
    currTime = '';
end

for j = 1:3 % grid search-dhs-ip
    if j == 1
        method = 'gs';
    elseif j == 2
        method = 'dhs';
    elseif j == 3
        method = 'ip';
    end
    % Write to File
    % Open Log File and get identifier
    logFile = [stationName,'.', currTime, method,'.log']; % output: log file name
    fID = fopen(fullfile(logFileFolder, logFile), 'wt');
    for i = 1:3 % E-N-U
%         writeInputLog(fID, stationName, coordinateName{i}, data{:, 'date'}, ...
%             polynDeg(i), osc{i}, heavJumps{i}, ...
%             resStor{i}(j).model.tst, resStor{i}(j).model.tau, resStor{1}(1).model.tstype, ...
%             KK, p, outlFactor);
        writeOutputLog(fID, dataStation, coordinateName{i}, ...
            resStor{i}(j).optp.polyn, ...
            reorderSinCos( resStor{i}(j).optp.oscil ), ...
            resStor{i}(j).optp.jumps, ...
            resStor{i}(j).model.tst, ...
            resStor{i}(j).optp.tsamp, ...
            resStor{i}(j).model.tau, ... 
            resStor{i}(j).error{1,2}, resStor{i}(j).error{2,2}...
            )
    end
    fclose(fID); % close log file
    fprintf('Calculation finished.\nPlotting and writing results (%s)...\n', method)
    
    % use parameters in file name
    resultSaveFile = fullfile(logFileFolder, [stationName, ...
        sprintf('_itrf%d_KK%d_p%.1f_outl%d', doITRFjump,KK, p, outlFactor), ...
        '.', method, ...
        '.csv']); % output: file name of computed trends (csv)
    
    % use custom fct
%     resultM = writeResultMatrix(tInterpol, trenddata, doITRFjump, KK, p, outlFactor, ...
%         [optParams{1}(1), optParams{2}(1), optParams{3}(1)], ...
%         [optParams{1}(2), optParams{2}(2), optParams{3}(2)]);
    
    % write matrix to csv file
    %writematrix(resultM, resultSaveFile, 'Delimiter', 'comma') % R2019a
    if doSaveResults; csvwrite(resultSaveFile, resultM); end% R2006
    
    % Visualize Results
    % set up title
    titleString = cell(3, 1); % preallocate
    for i = 1:3 % E-N-U
        % set up plot title
        if isempty(tauTypes{i});tsStr='none'; else; tsStr=sprintf('%s %s', transientType{i,1},transientType{i,2}); end
        titleString{i} = sprintf('Station:"%s" Transient:%s jump(itrf):%s  jump(eq):%s  RMS=%.2fmm WRMS=%.2fmm (%s)', ...
            stationName, tsStr,...
            mat2str(doITRFjump(i)), ...
            mat2str(doEQjump(i)), ...
            resStor{i}(j).error{1,2}, resStor{i}(j).error{2,2}, ... % rms,wrms
            method); 
    end
    
    % visualize time series and results
    figTSA = figure;
    VisualizeTS_Trend_Outliers_ITRF_ENU(...
        data{:, 'date'}, [data{:, 3}, data{:, 4}, data{:, 5}], outlier_logicals, coordinateName, ...
        tInterpol, [resStor{1}(j).trend, resStor{2}(j).trend, resStor{3}(j).trend], ...
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
        data{:, 3}-resStor{1}(j).trend, ... % E residuals
        data{:, 4}-resStor{2}(j).trend, ... % N residuals
        data{:, 5}-resStor{3}(j).trend], ...% U residuals
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
end
%%
fprintf('Done!\n')
% close all

function res = getParamStruct(params, os, nPoly, nOsci, nJump, nTran)
%gets result from parameter vector and writes them to struct.
s(1) = os+1;
s(2) = s(1)+nPoly;
s(3) = s(2)+nOsci;
s(4) = s(3)+nJump;
s(5) = s(4)+nTran;

res.error =  params( 1    : s(1)-1 );
res.polyn =  params( s(1) : s(2)-1 );
res.oscil =  params( s(2) : s(3)-1 );
res.jumps =  params( s(3) : s(4)-1 );
res.tsamp =  params( s(4) : s(5)-1 );
end

function res = reorderSinCos(osc)
%reorders sine/cosine oscillation components
% input: oscillation components as vector
% (number of coefficients must be even)
if mod(osc,2) == 1
    error('reorderSinCos: input error: oscillation component vector must contain even number of elements')
end
res = [osc(1:2:end); osc(2:2:end)];
end

function out = getTrendError(x, b, polynDeg, W, j_t, ts_t, tau, tsType, KK, p, outl_factor)
[~,results,~,~] = computeTrendIRLS(x, b, polynDeg, W, j_t, ts_t, tau, tsType, KK, p, outl_factor);
out = results{cellfun(@(x) strcmp(x,'rms'), results(:,1)),2};
end