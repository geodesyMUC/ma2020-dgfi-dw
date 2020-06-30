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

pol_N = length(pol); % number of polynomial coeffiecients
w_N = length(w); % number of periodic coefficients
jump_N = length(jt); % number of shifts
ts_N = length(ts_t); % number of transients

if length(tau) ~= length(ts_t) || length(ts_t) ~= length(tsType)
    error('LS error:transient model: length of tau,tau datetime and type of tau vectors do not match')
end

yy = zeros(size(x, 2), pol_N + w_N + jump_N + ts_N);

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
for i = 1:ts_N
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