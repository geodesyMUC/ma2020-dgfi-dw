function Y = TimeFunction(x, pol, CS, w, jt, b, eqjt, a, T)
% Creates time series from estimated parameters
% x: timestamps
% pol: polynomial coeffiecients
% CS: cos/sin periodic coefficients
% w: vector containing periods (rad)
% shift
% b: jump terms
% eq jump times
% amplitude of transient
% relaxation time T in seconds

if size(x, 1) > size(x, 2) % fix row/col 
    x = x';
end

pol_N = length(pol); % number of polynomial coeffiecients
w_N = length(w); % number of periodic coefficients
jump_N = length(jt); % number of shifts
eq_N = length(eqjt)*length(T); % number of logarithmic transients (n_events * n_tau)

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
% T = 1; % T = 1y -> constant
for i = 1:2:eq_N
    dt = x - eqjt(i);
    dt(dt < 0) = 0; % Every observation BEFORE the eq event
    % compute
    yy(:, cnt   ) = a(i  ) * log( 1 + dt./T(1) ); % logarithmic transient 1
%     yy(:, cnt) = a(i) * (1 - exp(-dt ./ T)); % exponential
    yy(:, cnt+1 ) = a(i+1) * log( 1 + dt./T(2) ); % logarithmic transient 2
    % increment counter with number of Tau
    cnt = cnt + length(T);
end

% x(t) = term1 + term2 + ... + termCNT
Y = sum(yy, 2); % row sum -> sum up all terms to compute y
end