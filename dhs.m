clear variables;
close all;
% down hill simplex algorithm

load('temp/output.mat')
tau1 = years(days(1:10:200));
tau2 = years(days(250:30:730));
lLim = [min(tau1), min(tau2)];
uLim = [max(tau1), max(tau2)];
tol = 0.0001;
x0 = [ min(tau1) + 0.0*(max(tau1)-min(tau1)) , min(tau2) + 0.0*(max(tau2)-min(tau2)) ];

figure
pcolor(tau1,tau2,-resultGrid);
hold on
shading interp;
contour(tau1,tau2,-resultGrid,'LineColor','k');
xlim([min(tau1) , max(tau1)])
ylim([min(tau2) , max(tau2)])
colorbar
hold on

nIter = 1000;
res = zeros(nIter,1);
for i=1:nIter
    res(i) = neldermead(x0,tau1,tau2,resultGrid,lLim,uLim,tol,false);
end
hold off

figure
histogram(res)

figure
pcolor(tau1,tau2,-resultGrid);
hold on
shading interp;
contour(tau1,tau2,-resultGrid,'LineColor','k');
xlim([min(tau1) , max(tau1)])
ylim([min(tau2) , max(tau2)])
colorbar
hold on
% ONeill
[xMin,fxMin,nFcall,nRest,err] = nelmin(...
    @myInterp2,...
    2, ...
    x0, ...
    tol, ...
    [max(tau1)-min(tau1), max(tau2)-min(tau2)], ...
    2, ...
    5000);

function out = neldermead(xInit,v1,v2,z,low,upp,tol,doPlot)
N = 2;
chi = 2;        % expansion
gamma = 0.5;    % contraction
sigma = 0.5;    % shrinkage/compression
rho = 1;        % reflection

p = 1/( N*sqrt(2) )*( N-1+sqrt( N+1 ) );
q = 1/( N*sqrt(2) )*( sqrt( N+1 )-1 );
L = [max(v1)-min(v1), max(v2)-min(v2)];

nF = 0;

if doPlot; plot(xInit(1),xInit(2),'kx'); end

x_ = zeros(3,2);
for i = 2:3 % vertices
    for j = 1:2 % dimensions
        if j == i-1
            x_(i,j) = xInit(j) + L(j)*p;
        elseif j ~= i-1
            x_(i,j) = xInit(j) + L(j)*q;
        end
%         % random scale factor
%         x_(i,j) = low(j) + rand(1,1)*(upp(j)-low(j));
    end
end
x_(1,:)=[];
X = [xInit; x_];

if doPlot
    pX = plotSimplex(X);
    pause(1)
    delete(pX)
end

% main iteration
s = Inf;
tries = 100;

while s > tol && tries >= 0
    fprintf('curr stddev: %.5f\n', s)
    if doPlot; pS = plotSimplex(X); end
    
    fX = interp2(v1,v2,z,X(:,1),X(:,2));
    nF = nF+1;
    [~,order] = sort(fX); % sort f(x) ascending
    X = X(order,:);
    
    % compute centroid of all vertices except x(end)
    m = sum( X(1:end-1,:) ,1) .* 1/N;
    m = shiftInBox(m,low,upp);
    if doPlot; pM = plot(m(1),m(2),'go'); end
    
    % reflect
    r = ( 1+rho ).*m - rho.*X(end,:);
    r = projectOnBounds(X(end,:), r, m, rho, low, upp);% boundary check
    fr = interp2(v1,v2,z,r(1),r(2));
    nF = nF+1;
    if doPlot; pR = plot(r(1),r(2),'mo'); end
    
    
    if fr < fX(1)
        % expand
        e = ( 1 + rho*chi ).*m - rho*chi.*X(end,:);
        e = projectOnBounds(X(end,:), e, m, rho*chi, low, upp); % boundary check
        fe = interp2(v1,v2,z,e(1),e(2));
        nF = nF+1;
        if doPlot; pE = plot(e(1),e(2),'cx'); end
        
        
        if fe < fr
            % accept e
            X(end,:) = e;
            fprintf('#%d:EXPANSION ACCEPTED\n', 101-tries);
        else
            % accept r
            X(end,:) = r;
            fprintf('#%d:REFLECTION ACCEPTED\n', 101-tries);
        end
        
    elseif fX(1) <= fr && fr < fX(end-1)
        % accept r
        X(end,:) = r;
        fprintf('#%d:REFLECTION ACCEPTED\n', 101-tries);
        
    elseif fX(end-1) <= fr && fr < fX(end)
        % outside contraction
        c = ( 1+ gamma*rho ).*m - gamma*rho*X(end-1,:);
        fc = interp2(v1,v2,z,c(1),c(2));
        nF = nF+1;
        if doPlot; pC = plot(c(1),c(2),'yx'); end
        
        if fc < fr
            % accept c
            X(end,:) = c;
            fprintf('#%d:OUTSIDE CONTRACTION ACCEPTED\n', 101-tries);
        else
            % shrink/compression
            X(2:end,:) = X(1,:) + sigma.*( X(2:end,:)-X(1,:) );
            fprintf('#%d:SHRINK', 101-tries);
        end
    else
        % inside contraction
        c = ( 1-gamma ).*m + gamma.*X(end,:);
        fc = interp2(v1,v2,z,c(1),c(2));
        nF = nF+1;
        if doPlot; pC = plot(c(1),c(2),'yx'); end
        if fc < fX(end)
            % if better than x(end), accept. if not, shrink
            X(end,:) = c;
            fprintf('#%d:INSIDE CONTRACTION ACCEPTED\n', 101-tries);
        else
            % shrink/compression
            X(2:end,:) = X(1,:) + sigma.*( X(2:end,:)-X(1,:) );
            fprintf('#%d:SHRINK\n', 101-tries);
        end
    end
    % calculate standard deviation of function values
    s = std(fX);
    % delete from plot
    if doPlot 
        pause(0.5)
        delete(pS)
        delete(pM)
        if exist('pC','var');delete(pC);end
        if exist('pE','var');delete(pE);end
        if exist('pR','var');delete(pR);end
    end
    tries = tries-1;
end
out = min(fX);
fprintf('result = %.4f\n', out);
end

function plotObj = plotSimplex(X)
    plotObj = plot( [X(:,1);X(1,1)] , [X(:,2);X(1,2)] ,'r');
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