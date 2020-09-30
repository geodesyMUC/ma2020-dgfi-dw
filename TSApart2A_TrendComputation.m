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

% inputFolder = 'station_data_dailyXYZfiles_xyz'; % XYZ files
inputFolder = 'station_data_dailyXYZfiles'; % Where Station Data (TSA_ReadAndTransform) is stored as ".mat"

jumpCSVLocation = 'src/jumps-ma.csv'; % Location of Jump Table/Jump Database
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
stationName = '21754S001A01'; % P-Okushiri - Hokkaido %[ok, 2 eqs, doeqjumps]
% stationName = '21778S001A01'; % P-Kushiro - Hokkaido %[ok, 2 eqs, doeqjumps]
% stationName = '23104M001A01'; % Medan (North Sumatra) %[ok, 2polynDeg, 2 eqs, doeqjumps]
% stationName = '41705M003A04'; % Santiago %[ok, doeqjumps]
% stationName = '41719M004A02'; % Concepcion %[ok]

% stationName = '41713S001A02'; % Los Angeles, Chile
% stationName = '21762S001A01'; % Kashiwazaki [3eq ok]
% stationName = '23114M001A01';  % Pulau Simuk, Ind

%%% Trend Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
polynDeg = [-1, -1, -1];% Polynomial Trend,Integer Degree, Range [-1..3]

% periods / oscillations in YEARS (=365.25 days) in vector form
osc = {[0.5 1], [0.5 1], [0.5 1]};     % common values: 0.5y, 1y
% osc = {[1], [1], [1]};     % common values: 0.5y, 1y
% osc = {[], [], []};

% Model ITRF jumps (set to "true") or ignore ITRF jumps (set to "false")
doITRFjump  = [false false false]; % E-N-U
doEQjump    = [true true true]; % E-N-U
doRemoveTs = false; % global
doTsOverlay = false; % global
estimationOpt = 1; % 1: E-N-U , 2: E&N-U , 3: E&N&U
tarFct = 'rms';

% specify type of transient: "log","exp","nil"
transientType = {...
    'log','exp'; ...    % coordinate1:E|X
    'log','exp'; ...    % coordinate2:N|Y
    'log','exp'};       % coordinate3:U|Z

% Parameter tau in [years] for computation of logarithmic transient for
% earthquake events (jumps):
% vector mapping different T (tau) relaxation coefficients [years]
tauVec1 = years(days(1:10:100));
tauVec2 = years(days(101:25:365*3));
% tauVec1 = years(days(1:5:365)); % full range tau

% optimization constraints for transient 1 and transient 2 [years]
lowLimit = [ 0.1/365.25, 20/365.25]; % simplex plot
uppLimit = [19/365.25, 5000/365.25]; % simplex plot

% lowLimit = [ 0.01/365.25, 250/365.25];
% uppLimit = [200/365.25, 8];

% weighting parameters
doWeighting = [false false false];
% control weight decay after eq.
wFactor = 1;  % [-1;1], -1 := strong decay ; 0 := linear decay ; 1 := no decay
twEq  = 0;      % time at which weighting takes eq weight "wEq" in [years], rel. time to ts
twNo = 2;       % time at which weighting takes default weight "wNo" in [years], rel. time to ts
wEq = 1;        % weight at time t_ts (=eq)
wNo = 0;        % default weight

% remove observations AFTER x years after eq [year]
removeAfter = 2;

% Additional Parameters for LSE/IRLSE
KK = 0;             % n of iterations for IRLS
p = 1.0;            % L_p Norm for IRLS
outlFactor = 5;     % median(error) + standard deviation * factor -> outlier

% other variables
coordinateName    = {'E [mm]', 'N [mm]', 'U [mm]'}; % used to label plots
% coordinateName    = {'X_{red} [mm]', 'Y_{red} [mm]', 'Z_{red} [mm]'}; % used to label plots
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
VisualizeTS_ENU2(data, dataStation, NaT(1,1) ...
    );
% VisualizeTS_ENU2(data, dataStation, currStationJumps{:, 2}, ...
%     'jumpTypes', currStationJumps{:, 4:5});

%% Prepare Data for LS Estimation
% Convert Oscillation Periods to Angular Velocity
oscW = cellfun(@(x) 2*pi./x, osc, 'UniformOutput', false);

% Time Series Timestamps
t  = data{:, 't'};       % timestamps [seconds];
t  = t ./  (86400 * 365.25) ; % convert timestamps to [JULIAN years]
% t  = years( seconds(t) ); % convert timestamps to [years]
t  = t-t(1);             % adjust for negative t values
t0 = data{1, 'date'};    % timestamp [datetime] of first measurement
tL = data{end, 'date'};  % timestamp [datetime] of last measurement

% Set up variables for relative date of jumps and transients
% Distinguish between EQs (invokes transient) and other jumps (do not invoke transient)
isEq = logical(currStationJumps{:, 4}); % 1:=earthquake; 0:=no earthquake

% Preallocate Storages
heavJumps  = cell(3,1); % Heaviside Jumps for the 3 coordinates E,N,U
transients = cell(3,1); % Transients for the 3 coordinates E,N,U

% Get Transient Datetime from jumps where type = eq
ts_t = getRelativeJumps_eq(currStationJumps{:, 2}, t0, tL, isEq);

% remove obs.
isValidIdx = zeros(size( t ));
ts_t2 = years( ts_t ) + removeAfter;
ts_t1 = years( ts_t );
for j = 1:length(ts_t)
    isValidIdx( t >= ts_t1(j) &  t < ts_t2(j) ) = 1;
end
t(~isValidIdx) = [];
t = t - t(1);
data(~isValidIdx, :) = [];
t0 = data{1, 'date'};   % Reassign
tL = data{end, 'date'}; % Reassign

ts_t = getRelativeJumps_eq(currStationJumps{:, 2}, t0, tL, isEq); % get updated transients (new t0, tL)

for i = 1:3 % E-N-U
    
    fprintf('[%s]: doEQjump set to "%s"\n', coordinateName{i}(1), mat2str(doEQjump(i)));
    if doEQjump(i)
        % all jump types
        jumps = getRelativeJumps(currStationJumps{:, 'Date'}, t0, tL);
    else
        % only jump types which are NOT eq
        jumps = getRelativeJumps_eq(currStationJumps{:, 2}, t0, tL, ~isEq);
    end
    
    % Add ITRF Realization Change Jumps
    fprintf('[%s]: doITRFjumps set to "%s"\n', coordinateName{i}(1), mat2str(doITRFjump(i)));
    if doITRFjump(i)
        % append itrf jumps to heaviside jumps vector
        jumps = [jumps; jumps0itrf];
    end
    
%     % Convert to [years] & assign to cell array 
%     transients{i} = years( ts_t );
%     heavJumps{i}  = years( jumps );
    
    % Convert to [JULIAN years] & assign to cell array 
    transients{i} = seconds( ts_t  ) ./  (86400 * 365.25) ;
    heavJumps{i}  = seconds( jumps ) ./  (86400 * 365.25);
    
    % Lookup Table for relative transient timestamps, types and limits (->dhs,ip)
    tsLUT{i} = getTransientReferences(transients{i}, transientType(i,:), ...
        lowLimit, uppLimit, doRemoveTs);
end

% create datetime array with equal date intervals (1d)
tInterpol =  data{:, 'date'}; % take original timestamps, no interpolation
% tInterpol =  min(data{:, 'date'}):days(1):max(data{:, 'date'}); % 1d interpolation
dateIntvlN = length(tInterpol); % n of timestamps

% preallocate output storages
resStor = cell(1,3); % E-N-U
nMethod = 3; % gs,dhs,ip
fields = {'optp', 'optx', 'model', 'trend', 'error', 'outl'};
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

%% Compute LS Weights

w = zeros(length(t), 3);
% figure % debug
for i = 1:3 % E,N,U
    if doWeighting(i)
        w(:,i) = wNo;
        tTs = unique( tsLUT{i}.time );
        
        figure % debug
        %         subplot(3,1,i) % debug
        for j = 1:length( tTs ) % Loop Transients
            twSta = [twEq; wEq];
            twMid = [twEq + (twNo-twEq)*0.5 ; wNo + (wEq-wNo)*0.5];
            
            if ( wEq-wNo ) < ( twNo-twEq )
                scale =  norm( twSta-twMid ) * ( wEq-wNo )/( twNo-twEq );
            else
                scale =  norm( twSta-twMid ) * ( twNo-twEq )/( wEq-wNo );
            end
            
            
            if wFactor >= 1
                warning('weighting: wFactor adjusted to be < 1')
                wFactor = 1-1e-5;
            elseif wFactor <= -1
                warning('weighting: wFactor adjusted to be > -1')
                wFactor = -1+1e-5;
            end
            % v1:
%             V = twSta - twMid;
%             nuV = [ V(2); -V(1) ] ./ norm( V ) * scale * wFactor;
%             twPivot = twMid + nuV;
            % v2:
            V = [twNo; wEq] - twMid;
            newV = V * wFactor;
            twPivot = twMid + newV;
            % ---
            dt = t - tTs(j);
            wLine = [...
                twEq, twPivot(1), twNo; ...
                wEq, twPivot(2), wNo];
            wx = interp1(wLine(1,:), wLine(2,:), dt, 'pchip', NaN);
            wxIdx = ~isnan(wx);
            w( wxIdx,i ) = wx( wxIdx );
            
            
            hold on % debug
            plot([twMid(1)+tTs(j) twPivot(1)+tTs(j)],[twMid(2) twPivot(2)], 'r'); % debug
            plot(twSta(1)+tTs(j), twSta(2), 'gx') % debug
            plot(twMid(1)+tTs(j), twMid(2), 'gx') % debug
            plot(twMid(1)+tTs(j), twMid(2), 'gx') % debug
        end
        plot(t,w(:,i),'b.') % debug
        
        grid on % debug
        ylabel('weight') % debug
        xlabel('t rel. to t_{0} [y]') % debug
        title(sprintf('Weighting for "%s" (%s)', coordinateName{i}(1), stationName)) % debug
        legend({'pivot line', 'weight controls'}, 'location', 'southoutside')
    else
        w(:,i) = deal(1);
    end   
end

%% Grid Search: generate tau vector
% for transients containing all combinations of values

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
        tauTypes{i} = {transientType{i,1}};
        tauVec      = tauVec1';
        tauCell{i}  = num2cell(tauVec,2);
        nTau(i)     = 1;
    elseif ~any(strcmp( transientType{i,1} , tsFctStr )) &&  any(strcmp( transientType{i,2} , tsFctStr ))
        % One Transient (tau_long)
        tauTypes{i} = {transientType{i,2}};
        tauVec      = tauVec2';
        tauCell{i}  = num2cell(tauVec,2);
        nTau(i)     = 1;
    else
        % No Transients, no match found
        tauVec = [];
        tauTypes{i} = {''};
        tauCell{i} = {[]};
        nTau(i) = 0;
    end
    % Lookup Table for relative transient timestamps, types and limits (->gs)
    % last argument needs to be false, transients must not be removed
    tsLUTgs{i} = getTransientReferences(transients{i}, tauTypes{i}, ...
        min(tauVec), max(tauVec), false);
end

% prepare array to store results rms,wrms,est.params FOR EVERY COMBINATION
% OF TAU (rows) and E,N,U (cols)
fxResAll = cell(3,1);
for i = 1:3 % E-N-U
    fxResAll{i} = zeros( size(tauCell{i},1) , 1 );
end

%% Grid Search: Set up Parameter Estimation Model
% other variables - get count of params (poly,osc, jumps, transients) per coordinate

for i = 1:3 % E-N-U
    nPolynTerms(i) = polynDeg(i) + 1; % 0, 1, 2, ...
    nOscParam(i)   = length(oscW{i}) * 2; % cos & sin components (C, S) for every oscillation
    nJumps(i)      = length(heavJumps{i}); % All Jumps - From DB and ITRF (if set to true)
    nEq(i)         = length(transients{i}); % n of eq events
    nTransients(i) = nEq(i) * nTau(i); % Only EQ Jumps -> n of transients
end

switch estimationOpt
    % outlier factor is first set to inf
    case 1 % E-N-U
        [params, computationName] = deal( cell(3,1) );
        for i = 1:3 % E-N-U
            params{i} = getParameterModel(KK, p, inf, t', data{:,2+i}', ...
                polynDeg(i), oscW{i}, heavJumps{i}', ...
                tsLUTgs{i}.time', tsLUTgs{i}.type', w(:,i)');
            computationName{i} = coordinateName{i};
        end
        
    case 2 % E&N-U
        if any( size(heavJumps{1}) ~= size(heavJumps{2}) )
            error('main: set up input parameter error: dim mismatch for heavJumps')
        elseif any( size(tsLUTgs{1}.time) ~= size(tsLUTgs{2}.time) )
            error('main: set up input parameter error: dim mismatch for ts times')
        elseif any( size(tsLUTgs{1}.type) ~= size(tsLUTgs{2}.type) )
            error('main: set up input parameter error: dim mismatch for ts types')
        end
        % E,N
        params{1} = getParameterModel(KK, p, inf, [t,t]', [data{:,3}, data{:,4}]', ...
            [polynDeg(1); polynDeg(2)], oscW{1}, [heavJumps{1}, heavJumps{2}]', ...
            [tsLUTgs{1}.time, tsLUTgs{2}.time]', [tsLUTgs{1}.type, tsLUTgs{2}.type]', [w(:,1), w(:,2)]');
        computationName{1} = [coordinateName{1}, ' & ', coordinateName{2}];
        
        % U
        params{2} = getParameterModel(KK, p, inf, t', data{:,5}', ...
            polynDeg(3), oscW{3}, heavJumps{3}', ...
            tsLUTgs{3}.time', tsLUTgs{3}.type', w(:,3)');
        computationName{2} = coordinateName{3};
    case 3 % E&N&U
        if any( size(heavJumps{1}) ~= size(heavJumps{2}) ) || any( size(heavJumps{2}) ~= size(heavJumps{3}) )
            error('main: set up input parameter error: dim mismatch for heavJumps')
        elseif any( size(tsLUTgs{1}.time) ~= size(tsLUTgs{2}.time) ) || any( size(tsLUTgs{2}.time) ~= size(tsLUTgs{3}.time) )
            error('main: set up input parameter error: dim mismatch for ts times')
        elseif any( size(tsLUTgs{1}.type) ~= size(tsLUTgs{2}.type) ) || any( size(tsLUTgs{2}.type) ~= size(tsLUTgs{3}.type) )
            error('main: set up input parameter error: dim mismatch for ts types')
        end
        % E,N,U
        params{1} = getParameterModel(KK, p, inf, [t,t,t]', [data{:,3}, data{:,4}, data{:,5}]', ...
            [polynDeg(1); polynDeg(2); polynDeg(3)], ...
            oscW{1}, ...
            [heavJumps{1}, heavJumps{2}, heavJumps{3}]', ...
            [tsLUTgs{1}.time, tsLUTgs{2}.time, tsLUTgs{3}.time]', ...
            [tsLUTgs{1}.type, tsLUTgs{2}.type, tsLUTgs{3}.type]', ...
            [w(:,1), w(:,2), w(:,3)]' );
        computationName{1} = [coordinateName{1}, ' & ', coordinateName{2},' & ', coordinateName{3}];
    otherwise 
        error('main: set up input parameter error: invalid value for estimationOpt')
        
end

%% Grid Search: Parameter Estimation

nComp = length(params); % n of computations
iRes = 1;               % count for processed coordinate components
for i = 1:nComp % E-N-U, E&N-U, E&N&U
    fxRes = zeros( size(tauCell{i},1) , 1 ); % preallocate
    
    % tau value combination loop (grid search)
    for j = 1:length(tauCell{i})
        % The following parameter tau need to be repeated for each transient or eq
        % to satisfy LS Input req.
        params{i}.tau = ...
            repmat( tauCell{i}{j} , [1,nEq(i)] ); % transient parameter tau
        
        [~, result_parameterC, ~, ~] = computeTrendIRLS( params{i}, doTsOverlay ); % LS
        fxRes(j) = result_parameterC{ strcmp(result_parameterC(:,1), tarFct) , 2};
        fxResAll{i}(j) = result_parameterC{ strcmp(result_parameterC(:,1), tarFct) , 2};
        
        fprintf('LS for "%s" (%d of %d), tau=[%s\b], rms=%.4f ...\n', ...
            computationName{i}, j, length(tauCell{i}), ...
            sprintf('%.4f ', tauCell{i}{j} ), fxRes(j));
    end
    [~, minIdx] = min( fxRes ); % minimization min( f(x) )
    minIdx = minIdx(1);         % prevent duplicates
    
    % set transient parameters to opt(x)
    params{i}.tau = repmat( tauCell{i}{minIdx} , [1,nEq(i)] );
    params{i}.outl = outlFactor; % assign outlier factor to remove outliers
    % LS: get parameters for opt(x)
    [~, result_parameterC, xEst, isOutlier] = computeTrendIRLS( params{i}, doTsOverlay );
    
    for j = 1:size(xEst,1)      % loop cells returned from LS
        if iscell(xEst)
            paramsEst = xEst{j}; % more than 1 component
        else
            paramsEst = xEst;    % exactly 1 component
        end
        optP = getParamStruct([result_parameterC{1,2}, result_parameterC{2,2}, paramsEst],...
            2, nPolynTerms(iRes), nOscParam(iRes), nJumps(iRes), nTransients(iRes));
        inputParams = getIdxParameterModel( params{i},j );
        
        optY = TimeFunction(... % calculate points
            inputParams.t, ...
            optP.polyn, ...
            optP.oscil, ...   % rearranged oscillation amplitudes (cosine, sine)
            inputParams.o, ...  % periods (global)
            inputParams.jt, ... % jump times
            optP.jumps, ...
            inputParams.tst, ...
            optP.tsamp, ...
            repmat( tauCell{iRes}{minIdx} , [1,nEq(iRes)] ), ... % opt tau (global)
            inputParams.tstype, ...
            doTsOverlay ...
            );
        
        resStor{iRes}(1).trend = optY;
        resStor{iRes}(1).error = result_parameterC;
        resStor{iRes}(1).optp  = optP;
        resStor{iRes}(1).optx  = tsLUTgs{iRes};
        resStor{iRes}(1).model = inputParams;
        resStor{iRes}(1).outl  = isOutlier(j,:);
        
        iRes = iRes+1; % increment coordinate component counter
    end
end

%% DHS,IP :Parameter Optimization Model

switch estimationOpt
    case 1 % E-N-U
        [lLim, uLim, x0, steps] = deal( cell(3,1) );
        for i = 1:3 % E-N-U
            params{i}.tst    = tsLUT{i}{:,'time'}';
            params{i}.tstype = tsLUT{i}{:,'type'}';
            params{i}.outl = inf;
            
            lLim{i}  = tsLUT{i}{:,'lBound'};
            uLim{i}  = tsLUT{i}{:,'uBound'};
            x0{i}    = tsLUT{i}{:,'lBound'}+1e-3;
            % steps for DHS init simplex,undefined for ts uBound==inf
            steps{i} = tsLUT{i}{:,'uBound'} - tsLUT{i}{:,'lBound'}; 
        end
    case 2 % E&N-U
        % E,N
        params{1}.tst    = [tsLUT{1}{:,'time'}, tsLUT{2}{:,'time'}]';
        params{1}.tstype = [tsLUT{1}{:,'type'}, tsLUT{2}{:,'type'}]';
        params{1}.outl = inf;
        % U
        params{2}.tst    = [ tsLUT{3}{:,'time'} ]';
        params{2}.tstype = [ tsLUT{3}{:,'type'} ]';
        params{2}.outl = inf;
        % E&N transients get assigned limits from first transient in LUT
        % -> they should be equal
        lLim{1}  = tsLUT{1}{:,'lBound'};
        uLim{1}  = tsLUT{1}{:,'uBound'};
        x0{1}    = tsLUT{1}{:,'lBound'}+1e-3;
        % steps for DHS init simplex,undefined for ts uBound==inf
        steps{1} = tsLUT{1}{:,'uBound'}-tsLUT{i}{:,'lBound'};
        
        lLim{2}  = tsLUT{3}{:,'lBound'};
        uLim{2}  = tsLUT{3}{:,'uBound'};
        x0{2}    = tsLUT{3}{:,'lBound'}+1e-3;
        % steps for DHS init simplex,undefined for ts uBound==inf
        steps{2} = tsLUT{3}{:,'uBound'}-tsLUT{i}{:,'lBound'};
        
    case 3 % E&N&U
        % E,N,U
        params{1}.tst    = [tsLUT{1}{:,'time'}, tsLUT{2}{:,'time'}, tsLUT{3}{:,'time'}]';
        params{1}.tstype = [tsLUT{1}{:,'type'}, tsLUT{2}{:,'type'}, tsLUT{3}{:,'type'}]';
        params{1}.outl = inf;
        % E&N&U transients get assigned limits from first transient in LUT
        % -> they should be equal
        lLim{1}  = tsLUT{1}{:,'lBound'};
        uLim{1}  = tsLUT{1}{:,'uBound'};
        x0{1}    = tsLUT{1}{:,'lBound'}+1e-3;
        % steps for DHS init simplex,undefined for ts uBound==inf
        steps{1} = tsLUT{1}{:,'uBound'}-tsLUT{i}{:,'lBound'};
        
    otherwise
        error('main: set up input parameter error: invalid value for estimationOpt (dhs,ip)')
end

%% DHS, IP: Optimization - Parameter Estimation
% & map plots from grid search results

for j = 2:3 % dhs-ip
    iRes = 1;
    for i = 1:nComp % E-N-U, E&N-U, E&N&U
        % GS Map Plots
        if  j == 2 && size(tauCell{i}{1},2) == 1% 1 tau
            figure
            plot(days(years( cell2mat(tauCell{i}) )) , fxResAll{i}(:) )
            xlabel('\tau [days]')
            ylabel('rms [mm]')
            title(["\tau parameter space for ", computationName{i}])
        elseif  j == 2 && size(tauCell{i}{1},2) == 2 && length(tauVec1)>1 && length(tauVec2)>1  % 2 tau
            resultGrid = reshape(fxResAll{i}(:,1), length( tauVec2 ) , length( tauVec1 ));
            
            if i == 1 % subplot. if not, just set up figure
                figure; %
                hold off % subplot
                subplot(1,3, i) % for subplots
            end % open only 1 figure at start

            figure
            
            xAxisPlot = days(years( tauVec1 )); % scale xaxis to days
            
            % v1 smooth gradient plot
            colormap(flipud(parula)) % low error = good
            pcolor( xAxisPlot , tauVec2, resultGrid);
            hold on
            shading interp;
            contour( xAxisPlot , tauVec2, resultGrid,'LineColor',[97, 97, 92]./255);
            
            % v2 distinct gradient plot
            % colormap(flipud(parula)) % low error = good
            % contourf(days(years(tauVec1)),days(years(tauVec2)),resultGrid);
            % contourf( days(years( tauVec1 )), tauVec2, resultGrid);
            
            % plot options
            xlabel('\tau_{1}short [days]')
            ylabel('\tau_{2}long [years]')
            title(['\tau parameter space for ', computationName{i}])
            c = colorbar;
            c.Label.String = 'RMS[mm]';
            xlim([min( xAxisPlot ) , max( xAxisPlot )])
            ylim([min(tauVec2) , max(tauVec2)])
            hold on
        end
                
        % Function for optimization
        optFun = @(x) getTrendError(...
            params{i}.t, ...     % t in years where t0=beginning of TS
            params{i}.b, ...     % vector with metric (coordinate)
            params{i}.w, ...     % weights
            params{i}.poly, ...  % polynome degree
            params{i}.o, ...     % periods
            params{i}.jt, ...    % jumps: time in years since t0
            params{i}.tst, ...   % eq transients: time in years since t0
            x, ...            % transient parameter tau
            params{i}.tstype,... % type (function) of tau 2(log|exp)
            params{i}.kk, ...    % n of iterations for IRLS
            params{i}.p, ...     % L_p Norm for IRLS
            params{i}.outl,...   % median(error) + standard deviation * factor -> outlier );
            doTsOverlay, ...  % overlay flag
            tarFct ...        % name of target fct for optimization
            );
        
        switch j
            case 2
                
%                 % DHS (custom method)
%                 restartScale = 0.001; %
%                 tol = 0.0001;
%                 doPrintDebug = true;
%                 [xMin, fxMin(i),nFnCalls(i),nRestarts(i),~] = ...
%                     dhscopt(optFun, x0{i}', steps{i}', lLim{i}', uLim{i}', tol, restartScale, doPrintDebug);

                % DHS (Bounded fminsearch)
                options = optimset('Display', 'final', 'TolX', 1e-6);
                xMin = fminsearchbnd(...
                    optFun, x0{i}, lLim{i}, uLim{i}, options...
                    );
                hold off
                
            case 3
                % MATLAB method
                if ~isempty(x0{i})
                    options = optimoptions(@fmincon,...
                        'Display','final','Algorithm','interior-point');
                    [xMin,fxMin_{i}] = fmincon(optFun,x0{i},...
                        [],[],[],[],lLim{i},uLim{i},[],options);
                else
                    xMin = [];
                end
        end
        % set transient parameters to opt(x)
        params{i}.tau = xMin;
        params{i}.outl = outlFactor; % assign outlier factor to remove outliers
        % LS: get parameters for opt(x)
        [~, result_parameterC, xEst, isOutlier] = computeTrendIRLS( params{i}, doTsOverlay );
        
        for k = 1:size(xEst,1)      % loop cells returned from LS
            if iscell(xEst)
                paramsEst = xEst{k};
            else
                paramsEst = xEst;
            end
            optP = getParamStruct([result_parameterC{1,2}, result_parameterC{2,2}, paramsEst],...
                2, nPolynTerms(iRes), nOscParam(iRes), nJumps(iRes), length(tsLUT{iRes}.time));
            inputParams = getIdxParameterModel( params{i},k );
            
            optY = TimeFunction(...     % calculate points:
                inputParams.t, ...
                optP.polyn, ...
                optP.oscil, ...         % rearranged oscillation amplitudes (cosine, sine)
                inputParams.o, ...      % periods (global)
                inputParams.jt, ...     % jump times
                optP.jumps, ...
                inputParams.tst, ...
                optP.tsamp, ...
                xMin, ...               % opt tau (global)
                inputParams.tstype, ...
                doTsOverlay ...
                );
            
            resStor{iRes}(j).trend = optY;
            resStor{iRes}(j).error = result_parameterC;
            resStor{iRes}(j).optp  = optP;
            resStor{iRes}(j).optx  = tsLUT{iRes};
            resStor{iRes}(j).model = inputParams;
            resStor{iRes}(j).outl  = isOutlier(k,:);
            iRes = iRes+1; % increment coordinate component counter
        end
    end
end

%% Evaluation
if ~doStaticFile
    currTime = [datestr(datetime('now'),'yyyy-mm-dd_HH-MM-SS'), '.'];
else
    currTime = '';
end

for j = 1:3 % grid search-dhs-ip
    switch j
        case 1
            method = 'gs';
        case 2
            method = 'dhs';
        case 3
            method = 'ip';
    end
    % Write to File
    % Open Log File and get identifier
    logFile = [stationName,'.', currTime, method,'.log']; % output: log file name
    
    fID = fopen(fullfile(logFileFolder, logFile), 'wt'); % log file ID
    writeHeaderLog(fID, stationName, data{1, 'date'}, data{end, 'date'}, length(data{:, 'date'}), ...
        estimationOpt, tarFct, doITRFjump, doTsOverlay, doRemoveTs, doWeighting, ...
        [wFactor, wEq, wNo, twEq, twNo])
    
    for i = 1:3 % E-N-U
        writeInputLog(fID, stationName, coordinateName{i}, data{:, 'date'}, ...
            polynDeg(i), osc{i}, heavJumps{i}, ...
            resStor{i}(j).model.tst, resStor{i}(j).optx.lBound, resStor{i}(j).optx.uBound, ...
            resStor{i}(j).model.tstype, ...
            KK, p, outlFactor);
        writeOutputLog(fID, dataStation, coordinateName{i}, ...
            resStor{i}(j).optp.polyn, ...
            resStor{i}(j).optp.oscil, ...
            resStor{i}(j).optp.jumps, ...
            resStor{i}(j).model.tst, ...
            resStor{i}(j).optp.tsamp, ...
            resStor{i}(j).model.tau, ...
            resStor{i}(j).outl, ...
            resStor{i}(j).error ...
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
        if isempty(tauTypes{i});tsStr='-'; else; tsStr=sprintf('%s-%s', transientType{i,1},transientType{i,2}); end
        
        % classic
%         titleString{i} = sprintf('Station:"%s" Transient:%s jump(itrf):%s  jump(eq):%s  RMS=%.2fmm WRMS=%.2fmm (%s)', ...
%             stationName, tsStr,...
%             mat2str(doITRFjump(i)), ...
%             mat2str(doEQjump(i)), ...
%             resStor{i}(j).error{1,2}, resStor{i}(j).error{2,2}, ... % rms,wrms
%             method);
        
        % thesis print
        titleString{i} =  sprintf('"%s" nP:%d nJ:%d nF:%d nEq:%d ts:%s RMS=%.2f WRMS=%.2f (%s)', ...
            stationName, ...
            polynDeg(i), length( jumps), length( osc{i} ), nEq(i), tsStr, ...   
            resStor{i}(j).error{1,2}, resStor{i}(j).error{2,2}, ... % rms,wrms
            method);
        
    end
    
    % visualize time series and results 
    figTSA = figure;
    VisualizeTS_Trend_Outliers_ITRF_ENU(...
        data{:, 'date'}, [data{:, 3}, data{:, 4}, data{:, 5}], ...
        [{resStor{1}(j).outl}; {resStor{2}(j).outl}; {resStor{3}(j).outl}]', ...
        coordinateName, ...
        tInterpol, [resStor{1}(j).trend, resStor{2}(j).trend, resStor{3}(j).trend], ...
        titleString, ...
        currStationJumps{:, 'Date'}, ...
        [currStationJumps{:, 'Earthquake'}, currStationJumps{:, 'HWSW_Change'}, currStationJumps{:, 'Unknown'}], ...
        jumpCategoryNames, ...
        []...%readITRFChanges(itrfChangesTextfile)...
        )
    % set(gcf, 'InnerPosition', [0 0 604 513]); % small figure
    set(gcf, 'InnerPosition', [0 0 1000 600]); % large figure
    
    % visualize residuals
    rsd(:,1) = data{:, 3}-resStor{1}(j).trend;
    rsd(:,2) = data{:, 4}-resStor{2}(j).trend;
    rsd(:,3) = data{:, 5}-resStor{3}(j).trend;
    
    figRes = figure;
    VisualizeResiduals(...
        data{:, 'date'}, ...
        rsd, ...% residuals
        [{resStor{1}(j).outl}; {resStor{2}(j).outl}; {resStor{3}(j).outl}]', ...
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
%% End
fprintf('Done!\n')
% close all

%% Functions
function res = getParamStruct(params, os, nPoly, nOsci, nJump, nTran)
%gets result from parameter vector and writes them to struct.
s(1) = os+1;
s(2) = s(1)+nPoly;
s(3) = s(2)+nOsci;
s(4) = s(3)+nJump;
s(5) = s(4)+nTran;

res.error = getParam(params, 1,    s(1)-1 );
res.polyn = getParam(params, s(1), s(2)-1 );
res.oscil = getParam(params, s(2), s(3)-1 );
res.jumps = getParam(params, s(3), s(4)-1 );
res.tsamp = getParam(params, s(4), s(5)-1 );

    function out = getParam(p, s, e)
        % p(arameters), s(tart index), e(nd index)
        if e>=s
            out = p(s:e);
        else
            out = [];
        end
    end
end

function res = getParameterModel(kk, p, outl, t, b, poly, o, jt, tst, tstype, w)
res.kk       = kk;      % n of iterations for IRLS
res.p        = p;       % L_p Norm for IRLS
res.outl     = outl;    % median(error) + standard deviation * factor -> outlier

% schema: rows are coordinates. 1 row = 1 coordinate, 2 rows = 2 c., ...
res.t        = fixRowCol( t );       % timestamps [years], must be sorted
res.b        = fixRowCol( b );       % vector with metric (coordinate)
res.poly     = poly;                 % integer polynome degree
res.o        = o;                    % periods [angular velocity]
res.jt       = jt;                   % jump timestamps [years] relative to t0
res.tst      = tst;                  % eq transients: time in years since t0
res.tstype   = tstype;               % type (function) of tau (log|exp)
res.w        = fixRowCol( w );       % weights

    function out = fixRowCol(in)
        if size(in,1) > size(in,2)
            out = in';
        else
            out = in;
        end
    end
end

function out = getIdxParameterModel(in, idx)
% INDEX: local=per coordinate, NO INDEX: global
out.kk       = in.kk;      % n of iterations for IRLS: NO INDEX
out.p        = in.p;       % L_p Norm for IRLS: NO INDEX
out.outl     = in.outl;    % median(error) + standard deviation * factor -> outlier: NO INDEX

% schema: rows are coordinates. 1 row = 1 coordinate, 2 rows = 2 c., ...
out.t        = in.t(idx, :);                    % timestamps [years], must be sorted: INDEX
out.b        = in.b(idx, :);                    % vector with metric (coordinate): INDEX
out.poly     = in.poly(idx, :);                 % integer polynome degree: INDEX
out.o        = in.o;                    % periods [angular velocity]: NO INDEX
out.jt       = in.jt(idx, :);                   % jump timestamps [years] relative to t0: INDEX
out.tst      = in.tst(idx, :);                  % eq transients: time in years since t0: INDEX
out.tstype   = in.tstype(idx, :);               % type (function) of tau (log|exp): INDEX
out.w        = in.w(idx, :);                    % weights: INDEX
if isfield( in, 'tau')
    out.tau      = in.tau;                          % tau: NO INDEX
end
end

function out = getTrendError(x, b, w, polynDeg, W, j_t, ts_t, tau, tsType, KK, p, outl_factor, doTsOverlay, tarFct)
[~,results,~,~] = computeTrendIRLS(x, b, w, polynDeg, W, j_t, ts_t, tau, tsType, KK, p, outl_factor, doTsOverlay);
if isempty(strcmp(x, tarFct))
    error('getTrendError: target function "%s" for optimization was not found', tarFct);
end
out = results{cellfun(@(x) strcmp(x,tarFct), results(:,1)),2};
end

