function [y, results, xEst, outlierLogical] = computeTrendIRLS(x, b, polynDeg, W, j_t, ts_t, tau, tsType, KK, p, outl_factor)
% IRLSE - Iterative Reweighted (Linear) Least Squares
% INPUT
%   x: vector containing time stamps in [SECONDS] relative to t0
%   b: vector containing observations
%   polynDeg: integer power of polynome denoting station velocity
%   w: vector containing periods in [RAD]
%   j_t: vector containing jump times in [SECONDS] relative to t0 (if not specified -> empty)
%   ts_t: vector containing transient times for earthquakes in [SECONDS] relative to t0 (if not specified -> empty)
%   (Note: t0 refers to the datetime of the first observation in the time series)
%   T: Logarithmic Transient Parameter in [YEARS], will be set to 1 if not specified
%   KK: number of iterations for IRLS
%   p: L_p Norm for IRLS. p=2 equals to the euclidean norm (=L2). 
%   If p=2, no reweighting will be applied, independent of the number of iterations KK
%   outl_factor: median(error)|mean(error) + standard deviation * factor -> outlier
doLog = false;

if nargin == 1
    % assume input is parameter struct
    b = x.b;                % 
    polynDeg = x.poly;      % 
    W = x.w;                % 
    j_t = x.jt;             % 
    ts_t = x.tst;           % 
    tau = x.tau;            % 
    tsType = x.tstype;      % 
    KK = x.kk;              % 
    p = x.p;                % 
    outl_factor = x.outl;   % 
    x = x.t;                % x needs to be reassigned
elseif nargin <= 8
    % Additional Parameters for (IR)LS, use default values
    KK = 0; % n of iterations for IRLS
    p = 2.0; % L_p Norm for IRLS
    outl_factor = 5; % median(error) + standard deviation * factor -> outlier
    if doLog;fprintf('No IRLS parameters defined. Calculating L2 norm LSE.\n');end
end

if doLog; fprintf('n of iterations for IRLS = %d,\np of L_p Norm for IRLS = %.2f,\nOutlier Factor = %d (mean of error + standard deviation * factor < outlier)\n', ...
    KK, p, outl_factor); end

if size(tau,2) ~= size(tsType,1)
    % try reshaping it (order of elements in tau is important!)
    % from [tau_1.1, tau_1.2, ... , tau_nEQ.1, tau_nEQ.2] to matrix
    if size(tau,2)/size(tsType,1) == length(ts_t) % assume ok: -> reshape
        % possible error: if they match by chance, computation will resume
        tau = reshape(tau, [length(ts_t), size(tsType,1)]);
        if size(tsType,1)>1
            tau = tau'; % so that rows->eq & cols->transients
        end
        doLocalTau = true;              % 1 set of tau for each eq (local)
    else % assume error: dim mismatch -> abandon
        error('transient error: n of tau per event does not match length of type of tau vector')
    end
elseif ~isempty(tau) % assume col count matches, continue row check 
    if size(tau,1) == 1                 % ok, 1 set of tau for all eq (global)
        doLocalTau = false;
    elseif size(tau,1) == length(ts_t)  % ok, 1 set of tau for each eq (local)
        doLocalTau = true;
    else                                % error
        error('transient error: n of tau does not match n of eq events OR is not 1')
    end
end

% convert datetimes from [seconds] to [years]
x = years(seconds(x));       % x./(365.25 * 86400);
j_t = years(seconds(j_t));      % jt./(365.25 * 86400);
ts_t = years(seconds(ts_t));    % eqjt./(365.25 * 86400);

% parameter counts
nData = length(x); % number of observations
nPolynTerms = polynDeg+1; % 0, 1, 2, ... 
nPeriodic = length(W); % oscillations
nPeriodicCoeff = nPeriodic*2; % cos & sin components (C, S) for every oscillation
nJumpCoeff = length(j_t); % All Jumps - From DB and ITRF (if set to true)
nEqParam = length(ts_t)*size(tau,2); % number of amplitudes for transients -> transient * number of eq

% Set up Map N: n of Parameter Storage Vector
N(1) = 0;
N(2) = N(1) + nPolynTerms;
N(3) = N(2) + nPeriodicCoeff;
N(4) = N(3) + nJumpCoeff;
N(5) = N(4) + nEqParam;

%% set up design matrix A
A = zeros(nData, N(5)); % initialize (measurements x unknown parameters)
if doLog; fprintf('Design Matrix A: %d x %d\n', size(A, 1), size(A, 2)); end

% 1:POLYNOMIAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 0:polynDeg
    A(:, N(1) + i + 1) = [x.^i]; % Needs +1 because of start at 0
end

% 2:OSCILLATION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 1; % counter for periodic coefficients

for i = 1:length(W)
    A(:, N(2)+cnt:N(2)+cnt+1) = [cos(x * W(i)), sin(x * W(i))];
    cnt = cnt+2; % increment to match coefficients
end

% 3:JUMP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nJumpCoeff
    % Heaviside Jump
    A(:, N(3) + i) = heaviside(x - j_t(i));
end

% 4:TRANSIENT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate logarithmic transient for all earthquakes in this TS
% For every EQ Event, there needs to be (n of Tau) columns
cnt=1;
for i = 1:length(ts_t)
    dt = x - ts_t(i);
    dt(dt < 0) = 0; % Every observation BEFORE the event
    if size(tau,2) > 0
        if doLocalTau; iEq = i; % local tau
        else; iEq = 1;          % global tau
        end
        
        % Choose function for Transient 1
        if strcmp(tsType(1,:),'log')
            A(:, N(4)+cnt  ) = log( 1 + dt./tau(iEq, 1) ); % logarithmic transient 1
        elseif strcmp(tsType(1,:),'exp')
            A(:, N(4)+cnt )  = exp(-dt./tau(iEq, 1));      % exponential transient 1
        end
        if size(tau,2) > 1
            % Choose function for Transient 2
            if strcmp(tsType(2,:),'log')
                A(:, N(4)+cnt+1 ) = log( 1 + dt./tau(iEq, 2) ); % logarithmic transient 2
            elseif strcmp(tsType(2,:),'exp')
                A(:, N(4)+cnt+1 ) = exp(-dt./tau(iEq, 2));      % exponential transient 2
            end
        end
    end
    cnt=cnt+size(tau,2);
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
    
    if doLog;fprintf('IRLS Iteration #%d\n', k);end
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

% wrap rms and wrms into result cell for output
results{1, 1} = 'rms';
results{1, 2} = rms;
results{2, 1} = 'wrms';
results{2, 2} = wrms;

if doLog;fprintf('WMRS = %.4f, RMS = %.4f\n', wrms, rms);end

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
y = TimeFunction(xSim, polynParam, periodicParam, W, j_t, jumpParam, ...
    ts_t, EQtransient, tau, tsType);
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
    if doLog; fprintf('Using Moore-Penrose pseudoinverse for Inversion of normal equation matrix.\n'); end
    x = pinv(A, 0.001) * B;
end
end