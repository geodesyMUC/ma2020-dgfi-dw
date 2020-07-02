function Y = TimeFunction(x, pol, CS, w, jt, b, ts_t, a, tau, tsType)
% Creates time series from estimated parameters
% x: timestamps [years]
% pol: polynomial coeffiecients
% CS: cos/sin periodic coefficients
% w: vector containing periods (rad)
% shift
% b: jump terms
% eq jump times [years]
% amplitude of transient
% relaxation time T [years]
% type of tau ("log"|"exp")

if size(x, 1) > size(x, 2) % fix row/col 
    x = x';
end

nPoly = length(pol); % number of polynomial coeffiecients
nW = length(w); % number of periodic coefficients
nJump = length(jt); % number of shifts
nTs = length(ts_t); % number of transients

if length(tau) ~= nTs || nTs ~= length(tsType)
    error('Time Function:transient model mismatch: vector length of tau,tau datetime and type of tau')
end

yy = zeros(size(x, 2), nPoly + nW + nJump + nTs);

cnt = 1; % counter variable for incremenation
% polynom terms
for i = 0:nPoly - 1
    yy(:, cnt) = pol(i + 1) * x.^(i);
    cnt = cnt + 1;
end

% periodic terms
for i = 1:nW
    yy(:, cnt) = CS(1, i) * cos(x * w(i)) + CS(2, i) * sin(x * w(i));
    cnt = cnt + 1;
end

% jump terms
for i = 1:nJump
    yy(:, cnt) = b(i) * heaviside(x - jt(i));
    cnt = cnt + 1;
end

% transient terms
for i = 1:nTs
    dt = x - ts_t(i); % num coeff!=num eq events
    dt(dt < 0) = 0; % Every observation BEFORE the eq event
    if strcmp(tsType(i),'log')
        yy(:,cnt)  = a(i) * log( 1 + dt./tau(i) ); % logarithmic transient
    elseif strcmp(tsType(i),'exp')
        yy(:,cnt) = a(i) * exp( -dt./tau(i) ); % exponential transient
    end
    cnt = cnt+1;
end

% x(t) = term1 + term2 + ... + termCNT
Y = sum(yy, 2); % row sum -> sum up all terms to compute y
end