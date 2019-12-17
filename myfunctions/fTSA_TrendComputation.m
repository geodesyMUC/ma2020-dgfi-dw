function [x] = fTSA_TrendComputation(stationname, jumpCSVPath, inputFolder, polynDeg, P, T, doITRFjump, imgDir)
% INPUT: Station Name (ie "AREQ")
% OUTPUT: Parameter Struct

% trend computation dgfi 17.6.2019 dw
% 5.8.2019 completed transient implementation
% 5.8.2019 rewriting as function with output 
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for station observation data (needs "<StationName>.mat"-files with 
% table containing 4 columns named "t", "E", "N", "U" and associated data)

% coordinate name vector
coordinateSTR = {'E', 'N', 'U'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select Station Data
% look in files
% pattern for data has to be "<StationName>.mat"
fpathpattern = fullfile(inputFolder, sprintf('%s.mat', stationname));
if exist(fpathpattern, 'file')
    load(fpathpattern)
else
    error('Specified station name .mat-file could not be found!')
end
% Reassign currStation variable
data = currStation.Data;

% Load Jump Table using custom function
dataJump = importfileJumpCSV(jumpCSVPath);
% Select only current station
TFJump = strcmp(dataJump{:, 'Station'}, stationname); % selection logical

% All Jumps for current station 
currStationJumps = dataJump(TFJump, :);
% All Jumps for current station only where Use is set to "1"
currStationJumps = currStationJumps(currStationJumps.Use == 1, :);

%% PREPARE PARAMETERS FOR TREND ESTIMATION
% convert time to YEARS (=365.25 days)
% reason: Numerical Stability in LSE matrix inversion
t = data{:, 't'} ./(365.25 * 86400);

% convert oscillations to angular velocity
W = 2 * pi ./ P;

% Jump table
% Convert jump to years: 0 = beginning of TS
t0 = data{1, 'date'};

% subtract time first TS observation to get relative times
jumps0 = etime(datevec(currStationJumps{:, 2}), ...
    datevec(repmat(t0, size(currStationJumps, 1), 1)));
% remove negative values -> jumps and eqs BEFORE Time Series starts
jumps0x = jumps0(jumps0 >= 0);
% convert jump times to years (365.25!!!! days) and sort asc
jumps0x = sort(jumps0x ./(365.25 * 86400));
HJump = [jumps0x];

% Distinguish between EQs (invokes log. transient) and other jumps (unknown, HW
% change -> do not invoke transient)
EQLogical = logical(currStationJumps{(jumps0 >= 0), 4}); % 1=eq, 0=other
EQHJump = [jumps0x(EQLogical)]; % EarthquakeHeavisideJump -> Calc Transient

%% Add ITRF Realization Change Jumps

ditrf(1) = datetime('2011-04-17', 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
ditrf(2) = datetime('2017-01-29', 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
% format yyyy-MM-dd
% extend here ...

if doITRFjump
    
    jumps0itrf = etime(datevec(ditrf), ...
        datevec(repmat(t0, size(ditrf, 1), 1)));
    % remove negative values -> itrf jumps BEFORE Time Series starts
    jumps0itrf = jumps0itrf(jumps0itrf >= 0);
    % convert jump times to years (365.25!!!! days) and sort asc
    jumps0itrf = sort(jumps0itrf ./(365.25 * 86400));
    
    % append to heaviside jumps vector
    HJump = [jumps0itrf; HJump]; % HeavisideJump
end

%% Prepare Trend Estimation
% variables
nPolynTerms = polynDeg + 1;
nOscParam = length(W) * 2;
nJumps = length(HJump); % All Jumps
nEQJumps = length(EQHJump); % Only EQ Jumps (Heaviside+Transient)

% preallocate (mostly needed for printing TS)
OUTLIERLOGICAL = cell(3, 1);
% equal intervals
% % dateIntvl =  data{:, 'date'}; % verify integrity of algorithm COMMENT
dateIntvl =  data{1, 'date'}:days(1):data{end, 'date'};
% n of intervals
dateIntvlN = length(dateIntvl);
trenddata = zeros(dateIntvlN, 3);

%% Trend Estimation with EQ Transients
for i = 1:3 % For E N U separately
    [y, rmse, xEst, outlierLogical] = computeTrendIRLS(t, data{:, i + 2}, ...
        polynDeg, ...  % polynome degree
        W, ...  % periods
        HJump, ...  % jump times
        EQHJump, ... % eq jumps
        T ... % transient parameter (default value 1y ok)
        );
    
    % split up xEst and add meta information
    x.coordinate(i, :) = coordinateSTR{i}; % E,N,U
    x.t0(i) = t0;
    x.ti(i) = data{end, 'date'};
    
    x.p(i, :) = xEst(1:nPolynTerms); % polynome coefficients
    
    x.osc_T(i, :) = P; % in [years]
    
    x.C(i, :) = xEst(nPolynTerms + 1:2:nPolynTerms + nOscParam); % C (oscillation)
    x.S(i, :) = xEst(nPolynTerms + 2:2:nPolynTerms + nOscParam); % S (oscillation)
    x.A(i, :) = sqrt(xEst(nPolynTerms + 1:2:nPolynTerms + nOscParam)'.^2 + ...
        xEst(nPolynTerms + 2:2:nPolynTerms + nOscParam)'.^2); % Amplitude A (oscillation)
    
    % maybe rework, temporary fix/workaround
    % happens if there are 
    if ~isempty(HJump)
        x.t_heaviside(i, :) = t0 + seconds(HJump * (365.25 * 86400));
        x.heaviside(i, :) = xEst(nPolynTerms + nOscParam + 1:nPolynTerms + nOscParam + nJumps);
    else
        x.t_heaviside(i, :) = zeros(0, 1);
        x.heaviside(i, :) = zeros(0, 1);
    end
    
    % maybe rework, temporary fix/workaround.
    % happens at station "BRAZ-54" among others. jumps but no EQ
    if ~isempty(EQHJump)
        x.t_logtrans(i, :) = t0 + seconds(EQHJump * (365.25 * 86400));
        x.logtrans(i, :) = xEst(...
            nPolynTerms + nOscParam + nJumps + 1:...
            nPolynTerms + nOscParam + nJumps + nEQJumps);
    else
        x.t_logtrans(i, :) = zeros(0, 1);
        x.logtrans(i, :) = zeros(0, 1);
    end
    OUTLIERLOGICAL{i} = outlierLogical; % logical containing outlier locations
    x.rmse(i, :) = rmse; % root square mean error
    trenddata(:, i) = y; % store computed trend points

end

% visualize function with EQs
% set up title
titleString = cell(3, 1);
titleStringPattern = 'Station "%s", n of Obs.=%d, ITRF Jumps: %s, RMSE=%.2fmm';

for i = 1:3 % for E N U respectively
    % Set up string based on doITRFjump
    if doITRFjump, ITRFstring = 'true';
    else, ITRFstring = 'false'; end
    
    titleString{i} = sprintf(titleStringPattern, stationname, size(data{:, 'date'}, 1), ...
        ITRFstring, x.rmse(i));
end

% Create invisible figure
fig = figure('visible','off');
% Call plot function
VisualizeTS_Trend_Outliers_ITRF_ENU(data{:, 'date'}, ...
    [data{:, 'E'}, data{:, 'N'}, data{:, 'U'}], ...
    dateIntvl, trenddata, ...
    titleString, ...
    OUTLIERLOGICAL, ...
    currStationJumps, ... % Jump table for this station
    ditrf); % ITRF datetime

set(gcf, 'InnerPosition', [0 0 1000 600]);

if ~exist(imgDir, 'dir')
    mkdir(imgDir) % Create Folder if it doesnt exist
    fprintf('Image Storage directory created.\n')
end

% Print Filenames & Save
fileNameImg = sprintf('%s-trend.png', stationname);
saveas(fig, fullfile(imgDir, fileNameImg));
close(fig);
end