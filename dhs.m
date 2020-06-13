clear variables;
close all;
% down hill simplex algorithm

load('temp/output.mat')
tau1 = years(days(1:10:200));
tau2 = years(days(250:30:730));
tol = 0.001;

figure
contourf(tau1,tau2,-resultGrid);
xlim([min(tau1) , max(tau1)])
ylim([min(tau2) , max(tau2)])
hold on

N = 2;
chi = 2;        % expansion
gamma = 0.5;    % contraction
sigma = 0.5;    % shrinkage/compression
rho = 1;        % reflection

% init simplex
% xInit = [ min(tau1) + (max(tau1)-min(tau1))/2 , min(tau2) + (max(tau2)-min(tau2))/2];
xInit = [ min(tau1) , min(tau2)];

p = 1/( N*sqrt(2) )*( N-1+sqrt( N+1 ) );
q = 1/( N*sqrt(2) )*( sqrt( N+1 )-1 );
L = [max(tau1)-min(tau1), max(tau2)-min(tau2)];

plot(xInit(1),xInit(2),'kx');

X = [xInit; zeros(3-1,2)];
for i = 2:3 % vertices
    for j = 1:2 % dimensions
        if j == i-1
            X(i,j) = xInit(j) + L(j);
        else
            X(i,j) = xInit(j);
        end
         
    end
end

pX = plotSimplex(X);
pause(1)
delete(pX)
% main iteration
done = false;
s = Inf;
tries = 100;

while s > tol && tries >= 0
    fprintf('curr stddev: %.5f\n', s)
    pS = plotSimplex(X);
    % sort f(x) ascending
    fX = interp2(tau1,tau2,resultGrid,X(:,1),X(:,2));
    [~,order] = sort(fX);
    X = X(order,:);
    
    % compute centroid of all vertices except x(end)
    m = sum( X(1:end-1,:) ,1) .* 1/N;
    pM = plot(m(1),m(2),'go');
    
    % reflect
    r = ( 1+rho ).*m - rho.*X(end,:);
    fr = interp2(tau1,tau2,resultGrid,r(1),r(2));
    pR = plot(r(1),r(2),'mo');
    
    if fr < fX(1)
        % expand
        e = ( 1 + rho*chi ).*m - rho*chi.*X(end,:);
        fe = interp2(tau1,tau2,resultGrid,e(1),e(2));
        pE = plot(e(1),e(2),'cx');
        
        if fe < fr
            % accept e
            X(end,:) = e;
        else
            % accept r
            X(end,:) = r;
        end
        
    elseif fX(1) <= fr && fr < fX(end-1)
        % accept r
        X(end,:) = r;
        
    elseif fX(end-1) <= fr && fr < fX(end)
        % outside contraction
        c = ( 1+ gamma*rho ).*m - gamma*rho*X(end-1,:);
        fc = interp2(tau1,tau2,resultGrid,c(1),c(2));
        pC = plot(c(1),c(2),'yx');
        
        if fc < fr
            % accept c
            X(end,:) = c;
        else
            % shrink/compression
            X(2:end,:) = delta.*X(1,:) + (1-delta).*X(2:end,:);
        end
    else
        % inside contraction
        c = ( 1-gamma ).*m + gamma.*X(end,:);
        fc = interp2(tau1,tau2,resultGrid,c(1),c(2));
        pC = plot(c(1),c(2),'yx');
        if fc < fX(end)
            % if better than x(end), accept. if not, shrink
            X(end,:) = c;
        else
            % shrink/compression
            X(2:end,:) = X(1,:) + sigma.*( X(2:end,:)-X(1,:) );
        end
    end
    % calculate standard deviation of function values
    s = std(fX);   
    % delete from plot
    pause(1)
    delete(pS)
    delete(pM)
    if exist('pC','var');delete(pC);end
    if exist('pE','var');delete(pE);end
    if exist('pR','var');delete(pR);end
    tries = tries-1;
end

function plotObj = plotSimplex(X)
    plotObj = plot( [X(:,1);X(1,1)] , [X(:,2);X(1,2)] ,'r');
end