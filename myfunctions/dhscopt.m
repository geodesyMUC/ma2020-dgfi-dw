function [xOpt, fxOpt, nF, restarts, err] = dhscopt(fn, xInit, L, low, upp, tol, restartScale)
err = '';

n = length(xInit);  % dimension of problem
chi = 2;            % expansion
gamma = 0.5;        % contraction
sigma = 0.5;        % shrinkage/compression
rho = 1;            % reflection

epsi = 0.001;       % default value for factorial test
% restartScale = 0.1; % default value for scale of restart simplex
nF = 0;             % number of fct calls

stepScale = 1;
p = 1/( n*sqrt(2) )*( n-1+sqrt( n+1 ) );    % initial simplex parameter 1
q = 1/( n*sqrt(2) )*( sqrt( n+1 )-1 );      % initial simplex parameter 2

maxRestarts = 5;
restarts = 0;
% Restart Loop
while restarts <= maxRestarts
    % Set up initial Simplex (axis by axis approach)
    x_ = [xInit; zeros(n+1-1, n)];
    for i = 2:n+1 % vertices
        for j = 1:n % dimensions
            if j == i-1
                x_(i,j) = xInit(j) + stepScale*L(j)*p;
            elseif j ~= i-1
                x_(i,j) = xInit(j) + stepScale*L(j)*q;
            end
            %         % random scale factor
            %         x_(i,j) = low(j) + rand(1,1)*(upp(j)-low(j));
        end
    end
    X = x_;

    % DHS Iteration Loop
    s = Inf;
    maxTries = 100;
    tries = maxTries;
    while s > tol && tries >= 0     
        fX = fn( X );
        nF = nF+1;
        [~,order] = sort(fX); % sort f(x) ascending
        X = X(order,:);       % adopt order
        
        % compute centroid of all vertices except x(end)
        m = sum( X(1:end-1,:) ,1) .* 1/n;
        m = shiftInBox(m,low,upp); % boundary check
        
        % reflect
        r = ( 1+rho ).*m - rho.*X(end,:);
        r = projectOnBounds(X(end,:), r, m, rho, low, upp); % boundary check
        fr = fn( r );
        nF = nF+1;
        
        if fr < fX(1)
            % expand
            e = ( 1 + rho*chi ).*m - rho*chi.*X(end,:);
            e = projectOnBounds(X(end,:), e, m, rho*chi, low, upp); % boundary check
            fe = fn( e );
            nF = nF+1;
            
            if fe < fr
                % accept e
                X(end,:) = e;
                % fprintf('#%d:EXPANSION ACCEPTED\n', 101-tries);
            else
                % accept r
                X(end,:) = r;
                % fprintf('#%d:REFLECTION ACCEPTED\n', 101-tries);
            end
            
        elseif fX(1) <= fr && fr < fX(end-1)
            % accept r
            X(end,:) = r;
            % fprintf('#%d:REFLECTION ACCEPTED\n', 101-tries);
            
        elseif fX(end-1) <= fr && fr < fX(end)
            % outside contraction
            c = ( 1+ gamma*rho ).*m - gamma*rho*X(end-1,:);
            fc = fn( c );
            nF = nF+1;
            
            if fc < fr
                % accept c
                X(end,:) = c;
                % fprintf('#%d:OUTSIDE CONTRACTION ACCEPTED\n', 101-tries);
            else
                % shrink/compression
                X(2:end,:) = X(1,:) + sigma.*( X(2:end,:)-X(1,:) );
                % fprintf('#%d:SHRINK', 101-tries);
            end
        else
            % inside contraction
            c = ( 1-gamma ).*m + gamma.*X(end,:);
            fc = fn( c );
            nF = nF+1;
            
            if fc < fX(end)
                % if better than x(end), accept. if not, shrink
                X(end,:) = c;
                % fprintf('#%d:INSIDE CONTRACTION ACCEPTED\n', 101-tries);
            else
                % shrink/compression
                X(2:end,:) = X(1,:) + sigma.*( X(2:end,:)-X(1,:) );
                % fprintf('#%d:SHRINK\n', 101-tries);
            end
        end
        % calculate standard deviation of function values
        s = std(fX);
        tries = tries-1;
    end
       
    % get minima for fx and x
    [fxOpt, fxOptIdx] = min(fX);
    xOpt = X(fxOptIdx,:);
    fprintf('result fx = %.4f\n', fxOpt);
    
    % RESTART if local minimum
    doRestart = false;
    
    for i = 1:n
        xTest = xOpt;
        delta = epsi * L(i); % step
        
        % Test 1 for local minimum
        xTest(i) = xTest(i) + delta;
        fxTest = fn(xTest);
        nF = nF+1;
        if fxTest < fxOpt
            doRestart = true;
            break
        end
        
        % Test 2 for local minimum
        xTest(i) = xTest(i) - delta - delta;
        fxTest = fn(xTest);
        nF = nF+1;
        if fxTest < fxOpt
            doRestart = true;
            break
        end
    end
    
    if doRestart
        % assume local minimum, restart
        % adapt stepScale, xInit
        xInit = xTest;
        stepScale = restartScale; % scale next initial simplex
        fprintf('RESTART\n');
        restarts = restarts + 1;
    else
        % assume global minimum, terminate
        break
    end
end

if restarts > maxRestarts
    err = sprintf('maximum number of restarts exceeded (%d). last best values returned', maxRestarts);
end

end

function p = projectOnBounds(x, p0, m, scale, lBound, uBound)
% random box method
% input vectors: investigated vertex, centroid, lower bounds, upper bounds
isLeft = p0 < lBound;
isRight = p0 > uBound;

% Following [Guin, 1968]
tries = 5;
if any(isLeft) || any(isRight)
    while tries >= 0 && (any(isLeft) || any(isRight))
        scale = scale/2;
        p = ( 1+scale ).*m - scale.*x;
        isLeft = p < lBound;
        isRight = p > uBound;
        tries = tries-1;
    end
    if tries == 0 % use random box method if Guin fails
        p = p0;
        
        a = 0.000001;
        b = 0.5;
        % project points outside lower boundary (RANDOM BOX)
        z = (b-a).*rand(1,nnz(isLeft)) + a;
        p(isLeft) = lBound(isLeft) + z.*( m(isLeft)-lBound(isLeft) );
        
        % project points outside upper boundary (RANDOM BOX)
        z = (b-a).*rand(nnz(isRight),1) + a;
        p(isRight) = uBound(isRight) - z.*( m(isRight)-lBound(isRight) );
    end
else
    p = p0;
end
end

function x = shiftInBox(x,lBound,uBound)
% box method
% input vectors: investigated vertex, lower bounds, upper bounds
x(x < lBound) = x(x < lBound)+0.000001;
x(x > uBound) = x(x > uBound)-0.000001;
end
