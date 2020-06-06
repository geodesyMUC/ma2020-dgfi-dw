function Y = TimeFunction(x, pol, CS, w, jt, b, ts, a, tau, tauT)
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
ts_N = length(ts)*length(tau); % number of transients (n_events * n_tau)

if length(tau) ~= size(tauT,1)
    error('length of tau vector for transients does not match length of type of tau vector')
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
% T = 1; % T = 1y -> constant
cnp=1; % counter for eq parameters
for i = 1:length(ts)%2:eq_N
    dt = x - ts(i); % num coeff!=num eq events
    dt(dt < 0) = 0; % Every observation BEFORE the eq event
    if ~isempty(tau)
        % compute
        if strcmp(tauT(1,:),'log')
            yy(:,cnt) = a(cnp) * log( 1 + dt./tau(1) ); % logarithmic transient 1
        elseif strcmp(tauT(1,:),'exp')
            yy(:, cnt)       = a(cnp  ) * exp( -dt./tau(1) ); % exponential 1
        end
        if length(tau)>1
            if strcmp(tauT(2,:),'log')
                yy(:, cnt+1 ) = a(cnp+1) * log( 1 + dt./tau(2) ); % logarithmic transient 2
            elseif strcmp(tauT(2,:),'exp')
                yy(:, cnt+1) = a(cnp+1) * exp( -dt./tau(2) ); % exponential 2
            end
        end
    end
    cnp = cnp + length(tau);
    cnt = cnt + length(tau); % increment counter with number of Tau
end

% x(t) = term1 + term2 + ... + termCNT
Y = sum(yy, 2); % row sum -> sum up all terms to compute y
end