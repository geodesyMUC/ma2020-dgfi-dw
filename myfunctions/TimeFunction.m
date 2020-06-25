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
ts_N = length(ts_t)*length(tau); % number of transients (n_events * n_tau)

if size(tau,2) ~= size(tsType,1)
    % try reshaping it (order of elements in tau is important!)
    % from [tau_1.1, tau_1.2, ... , tau_nEQ.1, tau_nEQ.2] to matrix
    if size(tau,2)/size(tsType,1) == length(ts_t) % assume ok: -> reshape
        % possible error: if they match by chance, computation will resume
        tau = reshape(tau, [length(ts_t), size(tsType,1)]);
        tau = tau'; % so that rows->eq & cols->transients
        doLocalTau = true;              % 1 set of tau for each eq (local)
    else % assume error: dim mismatch -> abandon
        error('transient error: n of tau per event does not match length of type of tau vector')
    end    
else % assume col count matches, continue row check 
    if size(tau,1) == 1                 % ok, 1 set of tau for all eq (global)
        doLocalTau = false;
    elseif size(tau,1) == length(ts_t)  % ok, 1 set of tau for each eq (local)
        doLocalTau = true;
    else                                % error
        error('transient error: n of tau does not match n of eq events OR is not 1')
    end
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
for i = 1:length(ts_t)
    dt = x - ts_t(i); % num coeff!=num eq events
    dt(dt < 0) = 0; % Every observation BEFORE the eq event
    if ~isempty(tau)
        if doLocalTau; iEq = i; % local tau
        else; iEq = 1;          % global tau
        end
        
        % compute
        if strcmp(tsType(1,:),'log')
            yy(:,cnt)  = a(cnp) * log( 1 + dt./tau(iEq,1) ); % logarithmic transient 1
        elseif strcmp(tsType(1,:),'exp')
            yy(:, cnt) = a(cnp) * exp( -dt./tau(iEq,1) ); % exponential 1
        end
        if size(tau,2)>1
            if strcmp(tsType(2,:),'log')
                yy(:, cnt+1) = a(cnp+1) * log( 1 + dt./tau(iEq,2) ); % logarithmic transient 2
            elseif strcmp(tsType(2,:),'exp')
                yy(:, cnt+1) = a(cnp+1) * exp( -dt./tau(iEq,2) ); % exponential 2
            end
        end
    end
    cnp = cnp + size(tau,2); % ???
    cnt = cnt + size(tau,2); % increment counter with number of Tau
end

% x(t) = term1 + term2 + ... + termCNT
Y = sum(yy, 2); % row sum -> sum up all terms to compute y
end