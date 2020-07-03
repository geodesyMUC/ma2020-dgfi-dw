function [y, results, xEst, outlierLogical] = computeTrendIRLS(x, b, polynDeg, W, j_t, ts_t, tau, tsType, KK, p, outl_factor, doTsOverlay)
% IRLSE - Iterative Reweighted (Linear) Least Squares
% INPUT
%   x: vector containing time stamps in [YEARS] relative to t0
%   b: vector containing observations
%   polynDeg: integer power of polynome denoting station velocity
%   w: vector containing periods in [RAD]
%   j_t: vector containing jump times in [YEARS] relative to t0 (if not specified -> empty)
%   ts_t: vector containing transient times for earthquakes in [YEARS] relative to t0 (if not specified -> empty)
%   (Note: t0 refers to the datetime of the first observation in the time series)
%   T: Logarithmic Transient Parameter in [YEARS], will be set to 1 if not specified
%   KK: number of iterations for IRLS
%   p: L_p Norm for IRLS. p=2 equals to the euclidean norm (=L2). 
%   If p=2, no reweighting will be applied, independently from the number of iterations KK
%   outl_factor: median(error)|mean(error) + standard deviation * factor -> outlier
doLog = false;

if nargin == 2
    % assume input is parameter struct: reassign
    doTsOverlay = b;        % overlay flag needs to be reassigned
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

if length(tau) ~= length(ts_t) || length(ts_t) ~= length(tsType)
    error('LS error:transient model: length of tau,tau datetime and type of tau vectors do not match')
end

A = createCoeffMat(x, polynDeg, W, j_t, ts_t, tau, tsType, doTsOverlay);

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

rms  = computeRMS(b - A*xEst);
wrms = computeWRMS(b - A*xEst, eye(length(b))); % weights 1 (diag matrix)

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
    
    rms = computeRMS(b - A*xEst);
    wrms = rms;%computeWRMS(b - A*xEst, w); % weights 1 (diag matrix)
    
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

y = A * xEst; % time series - trend
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

function [rms] = computeRMS(e)
nData = length(e); % number of observations
rms = sqrt(1/nData * sum((e).^2)); % root mean square
end

function [wrms] = computeWRMS(e, w_i)
wrms = sqrt( diag(w_i)'*e.^2 / sum( diag(w_i) ) ); % weighted root mean square
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