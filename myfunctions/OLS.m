function [y, results, xEst, outlierLogical] = OLS(x, b, w, polynDeg, osc, j_t, ts_t, tau, tsType, KK, p, outl_factor, doTsOverlay)
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
    doTsOverlay = b;        % overlay flag, needs to be reassigned
    b = x.b;                % 
    w = x.w;                %
    polynDeg = x.poly;      % 
    osc = x.o;              % 
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

% set up design matrix
nCoord = size(x,1);
nB = size(x,2);
A = zeros( nnz(x) , sum(polynDeg)+nCoord + length(osc)*2*nCoord + numel(j_t) + numel(ts_t) );
nxEst = 0;
nxEstStor = zeros( nCoord+1 , 1);

for i = 1:nCoord
    j = (i-1)*nB + 1;
    k = nxEst + 1;
    nxEst = nxEst + polynDeg(i)+1 + length(osc)*2 + size(j_t,2) + size(ts_t,2);
    nxEstStor(i+1) = nxEst;
    A( j:nB*i , k:nxEst ) = createCoeffMat(x(i,:), polynDeg(i), osc, j_t(i,:), ts_t(i,:), tau, tsType(i, :), doTsOverlay);
end

% timestamps, measurements, weights reshape
x = x';
x = x(:);

b = b';
b = b(:);

w = w';
w = w(:);

%% (1) Calculate initial parameters xEst from A, b (and weights w)
if w == ones(length(w),1)
    % OLS
    [xEst, e] = computeLeastSquares(A, b);
else
    % WLS
    [xEst, e] = computeWeightedLeastSquares(A, b, w);
%     xEst = lscov(A,b,w); % MATLAB method
    e = A*xEst - b;
end
%% (2) Detect & Remove outliers
% check if outliers are present

% v1
% outlierLogical = abs(e) > ( mean(e) + std(e)*outl_factor ); % Logical with outliers

% v2
% Logical with outliers based Median Absolute Deviation MAD
outlierLogical = abs(e) > ( abs( median(e) ) + median( abs(e - median(e)) ) * outl_factor); 

if nnz(outlierLogical) > 0
    fprintf('removed %s obs (MAD)\n', num2str(nnz(outlierLogical)));
end
% if so, then remove outliers and compute LSE one more time
if outl_factor < 100 && nnz(outlierLogical) > 0
    b(outlierLogical)    = []; % remove them from observation vector b
    A(outlierLogical, :) = []; % remove them from design matrix A
    w(outlierLogical)    = []; % remove them from weight vector w
    x(outlierLogical)    = []; % remove them from time vector t (x)
    % LS one more time
    [xEst, e] = computeLeastSquares(A, b);
end

%% compute errors
rms  = computeRMS(b - A*xEst);
wrms = computeWRMS(b - A*xEst, w); % weights 1 (diag matrix)

fnRms = @(idx) computeRMS(b(idx) - A(idx,:)*xEst);

rms35d = fnRms( getTimeIdx(years( days(35) ), x, ts_t) );
rms105d = fnRms( getTimeIdx(years( days(105) ), x, ts_t) );
rms1y = fnRms( getTimeIdx( 1, x, ts_t) );
rms2y = fnRms( getTimeIdx( 2, x, ts_t) );

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
e(e==0) = 1e-6; % add small value to errors == 0 to prevent div by 0
e = abs(e);

for k = 1:KK
    if doLog;fprintf('IRLS Iteration #%d\n', k);end
    w = e;
    
    [xEst, e] = computeWeightedLeastSquares(A, b, w);
    e(e==0) = 1e-6; % add small value to errors == 0 to prevent div by 0
    e = abs( e );
    
    rms = computeRMS(b - A*xEst);
    wrms = computeWRMS(b - A*xEst, e);
    
    e = abs(e).^(p - 2); % compute new weights
    
    %% Error each iteration (Used for Debugging/Monitoring iterative procedure)
    pNorm = norm(e, p); % error at each iteration
    E = [E pNorm];
    
    
    
    RMS_vector = [RMS_vector rms]; % append
    WRMS_vector = [WRMS_vector wrms]; % append
end

% %% (4) Debug Plots & Stats
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

% wrap errors into result cell for output
results{1, 1} = 'rms';
results{1, 2} = rms;
results{2, 1} = 'wrms';
results{2, 2} = wrms;
results{3, 1} = 'rms35d';
results{3, 2} = rms35d;
results{4, 1} = 'rms105d';
results{4, 2} = rms105d;
results{5, 1} = 'rms1y';
results{5, 2} = rms1y;
results{6, 1} = 'rms2y';
results{6, 2} = rms2y;

if doLog;fprintf('WMRS = %.4f, RMS = %.4f\n', wrms, rms);end

y = A * xEst; % time series - trend
xEst = xEst'; % return row
if nCoord > 1
    % preallocate
    [xEst_,  y_] = deal( cell(nCoord, 1) );
    % assign output to cells if nCoord is larger than 1
    outlierLogical = reshape( outlierLogical , [nB, nCoord] );
    outlierSum     = sum( ~outlierLogical ); % n of inliers
    outlierSum     = [0, outlierSum];
    for i = 1:nCoord
        outlierSum(i+1) = outlierSum(i)+outlierSum(i+1); % cumulative sum
        % part and reassign
        xEst_{i} = xEst( nxEstStor(i)+1 : nxEstStor(i+1) );
        y_{i}   = y( outlierSum(i)+1 : outlierSum(i+1) );
    end
    y = y_;
    xEst = xEst_;    
end
outlierLogical = outlierLogical'; % components -> rows, obs -> columns
end

%% Functions
function [xEst, e] = computeLeastSquares(A, b)
% Calculate parameters xEst and error vector e from A, b

% % Option 1 -----
% Nmat = A'*A; % normal equations
% n = A'*b;
% 

% Option 2 ----- Pseudoinverse (Penrose-Moore)
Nmat = A' * A; % normal equations
n = A' * b;
xEst = pinv(Nmat, 1e-8) * n;
% xEstx = Nmat\A'*b;

% % Option 3 --- Use Function
% xEst = lsqInvMMult(A' * A, A' * b);

% Detect & Remove outliers
e = A * xEst - b; % Error vector

end

function [xEst, e] = computeWeightedLeastSquares(A, b, w)
% Calculate parameters xEst and error vector e from A, b and the
% weight vector w

% w = abs(w).^(p - 2);       % compute Error weights w(i) [Wikipedia]
% denom_p = max(sum(w), 0.001); % Pick maximum to prevent zero denominator for new w(i)
% W = diag(w/denom_p);       % Normalized weight matrix

W = diag(w);

N = A'*W*A; % normal equations
if nnz(isnan(N))
    error('computeWeightedLeastSquares: error: matrix N = transp(A) * Wmat * A contains NaN.') % If matrix contains NaN -> Abort
end
n = A'*W*b;
% xEst = lsqInvMMult(N, n);
xEst = pinv(N, 1e-8) * n;

% WA = W * A; % [Burrows]
% N_ = WA' * WA; % [Burrows]

% xEst_ = pinv(N_) * (WA' * W) * b;
% xEst_ = lsqInvMMult(N_, (WA' * Wmat) * b); % [Burrows] % Alternative: Use Function

% xEst = pinv(N) * n;

e = A * xEst - b; % Error vector

end

function [rms] = computeRMS(e)
nData = length(e); % number of observations
rms = sqrt(1/nData * sum((e).^2)); % root mean square
end

function [wrms] = computeWRMS(e, w)
if size(w,1)>1 && size(w,2)>1
    w = diag(w); % make vector if matrix
end
wrms = sqrt( w'*e.^2 / sum( w ) ); % weighted root mean square
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
    x = pinv(A, 0.00001) * B;
end
end

function idx = getTimeIdx( tLim, t, ts )
%GETTIMEIDX creates logical vector
% input:
%   tLim time limit in [years]
%   t timestamp vector in [years]
%   ts transients in [years]
uniqt = unique( ts );
idx = false( length(t),1 );
for i = 1:length(uniqt)
    s = uniqt(i);
    e = uniqt(i) + tLim;
    idxNew = and( t>=s , t<e ) ;
    idx = or( idx, idxNew );
end
end