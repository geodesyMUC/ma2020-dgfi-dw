function [y, results, xEst, outlierLogical] = computeTrendIRLS(x, b, nPolyn, W, jt, eqjt, T, KK, p, outl_factor)
% IRLSE - 
% INPUT
%   x: vector containing time stamps in [SECONDS] relative to t0
%   b: vector containing observations
%   nPolyn: integer power of polynome denoting station velocity
%   w: vector containing periods in [RAD]
%   jt: vector containing jump times in [SECONDS] relative to t0 (if not specified -> empty)
%   eqjt: vector containing jump times for earthquakes in [SECONDS] relative to t0 (if not specified -> empty)
%   (Note: t0 refers to the datetime of the first observation in the time series)
%   T: Logarithmic Transient Parameter in [YEARS], will be set to 1 if not specified
%   KK: number of iterations for IRLS
%   p: L_p Norm for IRLS. p=2 equals to the euclidean norm (=L2). 
%   If p=2, no reweighting will be applied, independent of the number of iterations KK
%   outl_factor: median(error)|mean(error) + standard deviation * factor -> outlier

if nargin <= 7
    % Additional Parameters (default values)
    KK = 0; % n of iterations for IRLS
    p = 2.0; % L_p Norm for IRLS
    outl_factor = 4; % median(error) + standard deviation * factor -> outlier
    fprintf('No IRLS parameters defined. Calculating L2 norm LSE.\n')
end

fprintf('n of iterations for IRLS = %d,\np of L_p Norm for IRLS = %.2f,\nOutlier Factor = %d (mean of error + standard deviation * factor < outlier)\n', ...
    KK, p, outl_factor);

if nargin == 6 || nargin == 9
    T = 1; % assign default value 1y for logar. transient parameter T
end

% convert datetimes from [seconds] to [years]
x = years(seconds(x));       % x./(365.25 * 86400);
jt = years(seconds(jt));      % jt./(365.25 * 86400);
eqjt = years(seconds(eqjt));    % eqjt./(365.25 * 86400);

% parameter counts
nData = length(x); % number of observations
nPolynTerms = nPolyn+1; % 0, 1, 2, ... 
nPeriodic = length(W); % oscillations
nPeriodicCoeff = nPeriodic*2; % cos & sin components (C, S) for every oscillation
nJumpCoeff = length(jt); % All Jumps - From DB and ITRF (if set to true)
nEqParam = length(eqjt)*length(T); % number of eq jumps -> transient * number of T

% Set up Map N: n of Parameter Storage Vector
N(1) = 0;
N(2) = N(1) + nPolynTerms;
N(3) = N(2) + nPeriodicCoeff;
N(4) = N(3) + nJumpCoeff;
N(5) = N(4) + nEqParam;

%% set up design matrix A
A = zeros(nData, N(5)); % initialize (measurements x unknown parameters)
fprintf('Design Matrix A: %d x %d\n', size(A, 1), size(A, 2));

% 1:POLYNOMIAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 0:nPolyn
    A(:, N(1) + i + 1) = [x.^i]; % Needs +1 because of start at 0
end

% 2:OSCILLATION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0; % counter for periodic coefficients

for i = 1:length(W)
    A(:, N(2) +  1 + cnt:N(2) +  cnt + 2) = [cos(x * W(i)), sin(x * W(i))];
    cnt = cnt + 2; % increment to match coefficients
end

% 3:JUMP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nJumpCoeff
    % Heaviside Jump
    A(:, N(3) + i) = heaviside(x - jt(i));
end

% 4:TRANSIENT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate logarithmic transient for all earthquakes in this TS
% For every EQ Event, there needs to be (n of Tau) columns
for i = 1:2:nEqParam
    dt = x - eqjt(i);
    dt(dt < 0) = 0; % Every observation BEFORE the event
    % Transient 1
    A(:, N(4)+i  ) = log( 1 + dt./T(1) ); % logarithmic transient 1
%     A(:, N(4) + i) = 1 - exp(-dt ./ T(1)); % exponential transient 1
    % Transient 2
    A(:, N(4)+i+1) = log( 1 + dt./T(2) ); % logarithmic transient 2
end

%% (1) Calculate initial parameters xEst from A, b

[xEst, e] = computeLeastSquares(A, b);

%% (2) Detect & Remove outliers

% check if outliers are present
outlierLogical = abs(e) > mean(e) + std(e) * outl_factor; % Logical with outliers

% if so, then remove outliers and compute LSE one more time
if nnz(outlierLogical) > 0
    b(outlierLogical) = []; % remove them from observation vector b
    A(outlierLogical, :) = []; % remove them from design matrix A
    % LS one more time
    [xEst, e] = computeLeastSquares(A, b);
end  

rms = computeRMS(A, b, xEst);
wrms = computeWRMS(A, b, xEst, eye(length(b))); % weights 1 (diag matrix)

%% (3) IRLS - weight optimization
% to prevent unwanted behaviour when computing trends for TS with multiple
% jumps close to each other in time, apply the xEst parameters for jumps to
% trend, then compute new xEst

% WORK IN PROGRESS @27.1.2020 (Update 2/2020 - doesnt occur anymore)

% debugging/irls algorithm monitoring values %%%%%%%%%%%%%%%%%%%%%%
E = [];
RMS_vector = [rms];
WRMS_vector = [wrms];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:KK
    
    fprintf('IRLS Iteration #%d\n', k);
    e(e==0) = 0.0001; % add small value to 0, so div by 0 is prevented
    [xEst, e, w] = computeWeightedLeastSquares(A, b, e, p);
    
    %% Error each iteration (Used for Debugging/Monitoring iterative procedure)
    pNorm = norm(e, p); % error at each iteration
    E = [E pNorm];
    
    rms = computeRMS(A, b, xEst);
    wrms = computeWRMS(A, b, xEst, w); % weights 1 (diag matrix)
    
    RMS_vector = [RMS_vector rms]; % append
    WRMS_vector = [WRMS_vector wrms]; % append
    
end

% %% (4) QA Plots & Stats
% figure
% plot(0:KK, RMS_vector, 'bx-')
% hold on
% plot(0:KK, WRMS_vector, 'mx-')
% grid on
% title('rms/wrms error')
% xlabel('# Iteration ->')
% ylabel('error [mm]')
% 
% figure
% plot(1:KK, E, 'bx-')
% grid on
% title(sprintf('irls error: p norm (p=%.1f)', p))
% ylabel('p norm [mm]')
% xlabel('# Iteration ->')

% get rms and wrms into result cell for output
results{1, 1} = 'rms';
results{1, 2} = rms;
results{2, 1} = 'wrms';
results{2, 2} = wrms;

fprintf('WMRS = %.4f, RMS = %.4f\n', wrms, rms);

%% (5) sample equidistant values for TIME
% -> for time series with LSE estimated parameters

% Get Parameters and apply them on equally spaced time series

% Get Polynomial parameters
polynParam = xEst(N(1) + 1:N(2));
% Get periodic parameters (C&S: cos, sin)
periodicParam = xEst(N(2) + 1:N(3));
periodicParam = [periodicParam(1:2:end - 1)'; periodicParam(2:2:end)'];
% Get Jump/Unit Step/Heaviside Parameters
jumpParam = xEst(N(3) + 1:N(4));
% Get Jump/Unit Step/Heaviside Parameters
EQtransient = xEst(N(4) + 1:N(5));

% "simulated" time series (equally spaced measurements, depending on time interval)
tInterpolation = years(days(1)); % interpolation / sampling interval in YEARS
xSim = min(x) : tInterpolation : max(x); % interpolation
xSim = x'; % original values, no interpolation
% call function
% creates y values (e.g. up values) for "simulated" time series
y = TimeFunction(xSim, polynParam, periodicParam, W, jt, jumpParam, ...
    eqjt, EQtransient, T);
% y = A * xEst;
end

%% IRLS Custom Functions
function [xEst, e] = computeLeastSquares(A, b)
% Calculate parameters xEst and error vector e from A, b

% % Option 1 -----
% Nmat = A'*A; % normal equations
% n = A'*b;
% xEst = Nmat\A'*b;

% Option 2 ----- Pseudoinverse (Penrose-Moore)
Nmat = A' * A; % normal equations
n = A' * b;
xEst = pinv(Nmat) * n;

% % Option 3 --- Use Function
% xEst = lsqInvMMult(A' * A, A' * b);

% Detect & Remove outliers
e = A * xEst - b; % Error vector

end

function [xEst, e, w_i] = computeWeightedLeastSquares(A, b, e_t0, p)
% Calculate parameters xEst and error vector e from A, b, p and the
% previous error vector e_t0

w_i_ = abs(e_t0).^((p - 2)/2); % compute Error weights w(i) with /2 [Burrows]
w_i = abs(e_t0).^(p - 2); % compute Error weights w(i) [Wikipedia]
denom_p = max(sum(w_i), 0.001); % Pick maximum to prevent zero denominator for new w(i)
Wmat = diag(w_i/denom_p); % Normalized weight matrix

% % Using Moore Penrose Pseudoinverse with tolerance
% xEst3 = pinv(A' * Wmat * A, 0.001) * A' * Wmat * b; % Weighted LSE equation
% 
% % Using Moore Penrose Pseudoinverse
% xEst2 = pinv(A' * Wmat * A) * A' * Wmat * b; % Weighted LSE equation
% 
% % Using \
% xEst = A' * Wmat * A \ A' * Wmat * b; % Weighted LSE equation

Nmat = A' * Wmat * A; % normal equations
n = A' * Wmat * b;

WA = Wmat * A; % [Burrows]
N_ = WA' * WA; % [Burrows]
if nnz(isnan(N_))
    error('error: matrix N = transp(A) * Wmat * A contains NaN.') % Test: if matrix contains NaN -> Abort
end
xEst_ = pinv(N_) * (WA' * Wmat) * b;
% xEst_ = lsqInvMMult(N_, (WA' * Wmat) * b); % [Burrows] % Alternative: Use Function

xEst = pinv(Nmat) * n;
% xEst = lsqInvMMult(A' * A, A' * b); % Alternative: Use Function

e = A * xEst - b; % Error vector

end

function [rms] = computeRMS(A, b, xEst)
nData = length(b); % number of observations
rms = sqrt(1/nData * sum((b - A * xEst).^2)); % root mean square
end

function [wrms] = computeWRMS(A, b, xEst, w_i)
wrms = sqrt(sum(w_i .* (b - A * xEst).^2)/sum(w_i)); % weighted root mean square
end

function x = lsqInvMMult(A, B)
% Set warning as error temporarily, so try catch can be used
warning('error', 'MATLAB:nearlySingularMatrix');
warning('error', 'MATLAB:singularMatrix');

try
    x = A \ B;
catch
    % Using Moore Penrose Pseudoinverse with tolerance
    fprintf('Using Moore-Penrose pseudoinverse for Inversion of normal equation matrix.\n');
    x = pinv(A, 0.001) * B;
end
end