function [y, results, xEst] = computeNonlinearTrendLS(t, x0, l, nPolyn, W, jt, eqjt)
% convert datetimes to from [seconds] to [years]
t = years(seconds(t));       % x./(365.25 * 86400);
jt = years(seconds(jt));      % jt./(365.25 * 86400);
eqjt = years(seconds(eqjt));    % eqjt./(365.25 * 86400);

% parameter counts
nData = length(t); % number of observations
nPolynTerms = nPolyn + 1; % 0, 1, 2, ...
nPeriodic = length(W); % oscillations
nPeriodicCoeff = nPeriodic * 2; % cos & sin components (C, S) for every oscillation
nJumps = length(jt); % All Jumps - From DB and ITRF (if set to true)
nEqJumps = length(eqjt); % number of eq jumps -> transient

% Set up Parameter Storage Vector
N(1) = 0;
N(2) = N(1) + nPolynTerms;
N(3) = N(2) + nPeriodicCoeff;
N(4) = N(3) + nJumps;
N(5) = N(4) + nEqJumps;

%% set up design matrix A
A = zeros(nData, N(5)); % initialize (measurements x unknown parameters)
fprintf('Design Matrix A: %d x %d\n', size(A, 1), size(A, 2));

% 10 iterations at first
for i = 1:10
    % set up A
    
end
end