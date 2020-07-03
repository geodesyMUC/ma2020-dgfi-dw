function [A] = createCoeffMat(t, polyd, w, j_t, ts_t, tsTau, tsType, doTsOverlay)
%CREATECOEFFMAT Creates Coeffient Matrix A

%% prepare variables
% parameter counts
nData = length(t); % number of observations
nPolynTerms = polyd+1; % 0, 1, 2, ...
nPeriodic = length(w); % oscillations
nPeriodicCoeff = nPeriodic*2; % cos & sin components (C, S) for every oscillation
nJumpCoeff = length(j_t); % All Jumps - From DB and ITRF (if set to true)
nEqParam = length(ts_t); % number of amplitudes for transients -> transient * number of eq

% Set up Map N: n of Parameter Storage Vector
N(1) = 0;
N(2) = N(1) + nPolynTerms;
N(3) = N(2) + nPeriodicCoeff;
N(4) = N(3) + nJumpCoeff;
N(5) = N(4) + nEqParam;

%% set up design matrix
A = zeros(nData, N(5)); % initialize (measurements x unknown parameters)
% 1:POLYNOMIAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 0:polyd
    A(:, N(1) + i + 1) = [t.^i]; % Needs +1 because of start at 0
end

% 2:OSCILLATION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 1; % counter for periodic coefficients

for i = 1:length(w)
    A(:, N(2)+cnt:N(2)+cnt+1) = [cos(t * w(i)), sin(t * w(i))];
    cnt = cnt+2; % increment to match coefficients
end

% 3:JUMP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nJumpCoeff
    % Heaviside Jump
    A(:, N(3) + i) = heaviside(t - j_t(i));
end

% 4:TRANSIENT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate logarithmic transient for all earthquakes in this TS
for i = 1:length(ts_t)
    dt = t - ts_t(i);
    dt(dt < 0) = 0; % Every observation BEFORE the event
    % depending on flag, let transient overlay following earthquakes
    if ~doTsOverlay && i < length(ts_t) && ts_t(i+1)~=ts_t(i)
        dt(dt >= (ts_t(i+1)-ts_t(i)) ) = 0; % Every observation AFTER the NEXT event
    end
    
    if     strcmp(tsType(i),'log')
        A(:, N(4)+i  ) = log( 1 + dt./tsTau(i) ); % logarithmic transient
    elseif strcmp(tsType(i),'exp')
        A(:, N(4)+i )  = exp( -dt./tsTau(i) );    % exponential transient
    end
end
end

