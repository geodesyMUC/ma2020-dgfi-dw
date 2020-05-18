function [y, results, xEst] = computeNonlinearTrendLS(t, x0, b, L0, nPolyn, W, jt, eqjt, T)
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
nEqJumps = length(eqjt) * 2; % number of eq jumps -> transient AND tau_log

% Set up Parameter Storage Vector
N(1) = 0;
N(2) = N(1) + nPolynTerms;
N(3) = N(2) + nPeriodicCoeff;
N(4) = N(3) + nJumps;
N(5) = N(4) + nEqJumps;

% extend approximate values to adjust for tau_log approximates
eqCoeff = zeros(nEqJumps, 1);
cnt = 0; 
for i = 1:length(eqjt)
    eqCoeff(cnt + 1: cnt + 2) = [x0(N(4) + i), T]; % log coeff, tau_log approximate
    cnt = cnt + 2; % increment counter 
end
% reshape, update, extend approximate values vector x0
x0(N(4)+1 : N(4)+length(eqjt)) = [];
x0 = [x0; eqCoeff];

%% set up design matrix A
A = zeros(nData, N(5)); % initialize (measurements x unknown parameters)
fprintf('Design Matrix A: %d x %d\n', size(A, 1), size(A, 2));

% 10 iterations at first
for k = 1:50
%     fprintf('Iteration %i ...\n', k);
    
    % set up l
    l = b - L0;
    % set up A
    % 1:POLYNOMIAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 0:nPolyn
        A(:, N(1) + i + 1) = [t.^i]; % Needs +1 because of start at 0
    end
    
    % 2:OSCILLATION MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cnt = 0; % counter for periodic coefficients
    
    for i = 1:length(W)
        A(:, N(2) +  1 + cnt:N(2) +  cnt + 2) = [cos(t * W(i)), sin(t * W(i))];
        cnt = cnt + 2; % increment to match coefficients
    end
    
    % 3:JUMP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:nJumps
        % Heaviside Jump
        A(:, N(3) + i) = heaviside(t - jt(i));
    end
    
    % 4:TRANSIENT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate logarithmic transient and tau for all earthquakes in this TS
    cnt = 1; % counter for transient coefficients
    for i = 1:length(eqjt)
        dt = t - eqjt(i);
        dt( dt<0 ) = 0; % Every observation BEFORE the event
        a0 = x0( N(4)+cnt ); % approximate value for log amplitude coefficient
        tau0 = x0 ( N(4)+cnt+1 ); % approximate value for log relaxation coefficient
        % Transient
        A(: , N(4)+cnt : N(4)+cnt+1) = [...
            log(1 + dt./tau0), ...              % function derivate a
            -(a0*dt)./(tau0^2 + dt.*tau0) ...   % function derivate tau
            ];
        % A(:, N(4) + i) = 1 - exp(-dt ./ T); % exponential transient
        
        cnt = cnt + 2; % increment to match coefficients
    end
    
    [dxEst, e] = computeLeastSquares(A, l);
    if isreal(dxEst)    % check if complex numbers
        x0 = x0 + dxEst;
        L0 = A*x0;
    else                % if vector contains complex numbers: abandon and throw warning
        warning('NLLSE: estimated parameters contain complex numbers - iterating abandoned, current values returned.')
        break
    end
    
end

xEst = x0;
y = L0;

% compute results for evaluation
rms = computeRMS(A, b, xEst);
% wrms = computeWRMS(A, b, xEst, w); % weights 1 (diag matrix) % NO WRMS IMPLEMENTED YET
% get rms and wrms into result cell for output
results{1, 1} = 'rms';
results{1, 2} = rms;
results{2, 1} = 'wrms';
results{2, 2} = rms; % NO WRMS IMPLEMENTED YET

fprintf('WMRS = %.4f, RMS = %.4f\n', rms, rms);


end

%% NonlinearLS Custom Functions
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

function [rms] = computeRMS(A, b, xEst)
nData = length(b); % number of observations
rms = sqrt(1/nData * sum((b - A * xEst).^2)); % root mean square
end

function [wrms] = computeWRMS(A, b, xEst, w_i)
wrms = sqrt(sum(w_i .* (b - A * xEst).^2)/sum(w_i)); % weighted root mean square
end