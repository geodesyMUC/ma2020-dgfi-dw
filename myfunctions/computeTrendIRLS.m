function [y, results, xEst, outlierLogical] = computeTrendIRLS(x, b, nPolyn, w, jt, eqjt, T, KK, p, outl_factor)
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

fprintf('n of iterations for IRLS = %d,\np of L_p Norm for IRLS = %.2f,\nOutlier Factor = %d (median of error + standard deviation * factor < outlier)\n', ...
    KK, p, outl_factor);

if nargin == 6 || nargin == 9
    T = 1; % assign default value 1y for logar. transient parameter T
end

% convert datetimes to from [seconds] to [years]
x = x./(365.25 * 86400);
jt = jt./(365.25 * 86400);
eqjt = eqjt./(365.25 * 86400);

% parameter counts
nData = length(x); % number of measurements
nPolynTerms = nPolyn + 1; % 0, 1, 2, ... 
nPeriodic = length(w); % oscillations
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

% 1:POLYNOMIAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 0:nPolyn
    A(:, N(1) + i + 1) = [x.^i]; % Needs +1 because of start at 0
end

% 2:OSCILLATION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0; % counter for periodic coefficients

for i = 1:length(w)
    A(:, N(2) +  1 + cnt:N(2) +  cnt + 2) = [cos(x * w(i)), sin(x * w(i))];
    cnt = cnt + 2; % increment to match coefficients
end

% 3:JUMP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nJumps
    % Heaviside Jump
    A(:, N(3) + i) = heaviside(x - jt(i));
end

% 4:TRANSIENT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate logarithmic transient for all earthquakes in this TS
for i = 1:nEqJumps
    dt = x - eqjt(i);
    dt(dt < 0) = 0; % Every observation BEFORE the eq event
    % Transient
    A(:, N(4) + i) = log(1 + dt ./ T);
end

%% (1) Calculate initial parameters xEst from A, b

% % Option 1 -----
% Nmat = A'*A; % normal equations
% n = A'*b;
% xEst = Nmat\A'*b;

% % Option 2 ----- Pseudoinverse (Penrose-Moore)
% Nmat = A'*A; % normal equations
% n = A' * b;
% xEst = pinv(Nmat) * n;

% Option 3 --- Use Function
xEst = lsqInvMMult(A' * A, A' * b);

%% (2) Detect & Remove outliers
e = A * xEst - b; % Error vector

% check if outliers are present
outlierLogical = abs(e) > median(e) + std(e) * outl_factor; % Logical with outliers

% if so, then remove outliers and compute LSE one more time
if nnz(outlierLogical) > 0
    b(outlierLogical) = []; % remove them from measurements
    A(outlierLogical, :) = []; % remove them from design matrix
    
    % Option 1 -----
    % Nmat = A'*A; % normal equations
    % n = A'*b;
    % xEst = Nmat\A'*b;
    
    % % Option 2 ----- Pseudoinverse (Penrose-Moore)
    % xEst = pinv(Nmat)* A'* b;
    
    % Option 3 --- Use Function
    xEst = lsqInvMMult(A' * A, A' * b);
    
    % update error vector
    e = A * xEst - b;
    % update nData
    nData = length(b);
end  

%% (3) IRLS - weight optimization
% to prevent unwanted behaviour when computing trends for TS with multiple
% jumps close to each other in time, apply the xEst parameters for jumps to
% trend, then compute new xEst

% WORK IN PROGRESS @27.1.2020

% (IRLS: Iterative Reweighted Least Squares, C.Sidney Burrus)
w_i = ones(size(b, 1), 1);

% compute WEIGHTED RMS and RMS error
error_pnorm = norm(e, p);

% wrmse_ = sqrt(1/nData * sum(w_i .* (b - A * xEst).^2));
wrmse_ = sqrt(sum(w_i .* (b - A * xEst).^2)/sum(w_i));
rmse_ = sqrt(1/nData * sum((b - A * xEst).^2));

% debugging/irls algorithm monitoring values %%%%%%%%%%%%%%%%%%%%%%
E = [];
RMS_ = [rmse_];
WRMS_ = [wrmse_];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:KK
    
    fprintf('IRLS Iteration #%d\n', k);
    e = A * xEst - b; % Error vector
    e(e==0) = 0.0001; % add small value to 0, so div by 0 is prevented
    w_i = abs(e).^((p - 2)/2); % compute Error weights w(i) with /2
    denom_p = max(sum(w_i), 0.001); % Pick maximum to prevent zero denominator for new w(i)
    Wmat = diag(w_i/denom_p); % Normalized weight matrix
    
%     % Using Moore Penrose Pseudoinverse with tolerance
%     xEst3 = pinv(A' * Wmat * A, 0.001) * A' * Wmat * b; % Weighted LSE equation

%     % Using Moore Penrose Pseudoinverse
%     xEst2 = pinv(A' * Wmat * A) * A' * Wmat * b; % Weighted LSE equation

%     % Using \
%     xEst = A' * Wmat * A \ A' * Wmat * b; % Weighted LSE equation

    % use custom Function for MMult
    % if matrix is close to being singular (inversion possibly not
    % possible), the pseudoinverse is used to do the matrix multiplication
    WA = Wmat * A;
    N_ = WA' * WA;
%     N_ = A' * Wmat * A;
    
    if nnz(isnan(N_))
        error('error: matrix N = transp(A) * Wmat * A contains NaN') % Test: if matrix contains NaN -> Abort
    end
    
%     xEst_  = lsqInvMMult(A' * Wmat * A, A' * Wmat * b); % old
    xEst = lsqInvMMult(WA' * WA, (WA' * Wmat) * b);
    
    %% Error each iteration (Used for Debugging/Monitoring iterative procedure)
    ee = norm(e, p); % error at each iteration
    E = [E ee];
    
    wrms_ = sqrt(sum(w_i .* (b - A * xEst).^2)/sum(w_i));
    rms_ = sqrt(1/nData * sum((b - A * xEst).^2));
    
    RMS_ = [RMS_ rms_];
    WRMS_ = [WRMS_ wrms_];
    
end


%% (4) QA Plots & Stats
figure
plot(0:KK, RMS_, 'bx-')
hold on
plot(0:KK, WRMS_, 'mx-')
grid on
title('rms/wrms error')
xlabel('# Iteration ->')
ylabel('error [mm]')

figure
plot(1:KK, E, 'bx-')
grid on
title(sprintf('irls error: p norm (p=%.1f)', p))
ylabel('p norm [mm]')
xlabel('# Iteration ->')

% compute WEIGHTED RMS and RMS error -> store in result cell
wrmse = sqrt(sum(w_i .* (b - A * xEst).^2)/sum(w_i));
rmse = sqrt(1/nData * sum((b - A * xEst).^2));

results{1, 1} = 'rms';
results{1, 2} = rmse;
results{2, 1} = 'wrms';
results{2, 2} = wrmse;

fprintf('WMRS = %.4f, RMS = %.4f\n', wrmse, rmse);

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

% "simulated" time series (equally spaced, depending on time interval)
xSim = x(1):1/365.25:x(end); 

% call function
% creates y values (e.g. up values) for "simulated" time series
y = TimeFunction(xSim, polynParam, periodicParam, w, jt, jumpParam, ...
    eqjt, EQtransient);
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

function Y = TimeFunction(x, pol, CS, w, jt, b, eqjt, a)
% Creates time series from estimated parameters
% x: timestamps
% pol: polynomial coeffiecients
% CS: cos/sin periodic coefficients
% w: vector containing periods (rad)
% shift
% logarithmic transient

pol_N = length(pol); % number of polynomial coeffiecients
w_N = length(w); % number of periodic coefficients
jump_N = length(jt); % number of shifts
eq_N = length(eqjt); % number of logarithmic transients

yy = zeros(size(x, 2), pol_N + w_N + jump_N + eq_N);

cnt = 1; % counter variable for incremenation
% polynom terms
for i = 0:pol_N - 1
    yy(:, cnt) = pol(i + 1) * x.^(i);
    cnt = cnt + 1;
end

% periodic terms
for i = 1:w_N
    yy(:, cnt) = CS(1, i) * cos(x * w(i)) + CS(2, i) * sin(x * w(i));
    cnt = cnt + 1;
end

% jump terms
for i = 1:jump_N
    yy(:, cnt) = b(i) * heaviside(x - jt(i));
    cnt = cnt + 1;
end

% transient terms
T = 1; % T = 1y -> constant
for i = 1:eq_N
    dt = x - eqjt(i);
    dt(dt < 0) = 0; % Every observation BEFORE the eq event
    % compute
    yy(:, cnt) = a(i) * log(1 + dt ./ T);
    % increment counter
    cnt = cnt + 1;
end

% x(t) = term1 + term2 + ... + termCNT
Y = sum(yy, 2); % row sum -> sum up all terms to compute y
end