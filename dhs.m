clear variables;
close all;
% down hill simplex algorithm

load('temp/output.mat')
tau1 = years(days(1:10:200));
tau2 = years(days(250:30:730));

figure
contourf(tau1,tau2,-resultGrid);
hold on

N = 2;
alpha = 1;      % reflection
gamma = 2;      % expansion
beta = 0.5;     % contraction
delta = 0.5;    % shrink

% init simplex
xInit = [ min(tau1) + (max(tau1)-min(tau1))/2, min(tau2) + (max(tau2)-min(tau2))/2];
plot(xInit(1),xInit(2),'kx');

scalefactor = 0.05;

X = [xInit; [xInit(1) + scalefactor*1,xInit(2)]; [xInit(1),xInit(2) + scalefactor*1]];
plotSimplex(X)

% main iteration
done = false;
tries = 20;
while done == false && tries >= 0
    pS = plotSimplex(X);
    % 2: SORT F(X) ASC
    FX = interp2(tau1,tau2,resultGrid,X(:,1),X(:,2));
    [~,order] = sort(FX);
    X = X(order,:);
    
    % 3: CALC CENTROID
    m = sum( X(1:end-1,:) ,1) .* 1/N;
    pM = plot(m(1),m(2),'go');
    
    % 4: REFLECT
    r = (1+alpha).*m - alpha.*X(end,:);
    fr = interp2(tau1,tau2,resultGrid,r(1),r(2));
    pR = plot(r(1),r(2),'mo');
    
    % 5: DECIDE -> EXPAND?
    if fr < FX(1)
        % EXPAND
        e = (1+gamma).*m - gamma.*X(end,:);
        fe = interp2(tau1,tau2,resultGrid,e(1),e(2));
        pE = plot(e(1),e(2),'cx');
        % TAKE r OR e
        if fe < fr
            X(end,:) = e;
        else
            X(end,:) = r;
        end
    elseif fr < FX(end-1,:)
        % 6
        X(end,:) = r;
    else
        if fr < FX(end,:)
            X(end,:) = r;
            FX(end) = fr;
        end
        
        % 7: CONTRACTION
        c = beta.*m + (1-beta).*X(end,:);
        pC = plot(c(1),c(2),'yx');
        fc = interp2(tau1,tau2,resultGrid,c(1),c(2));
        
        % 8: DECIDE -> COMPRESS?
        if fc < FX(end,:)
            X(end,:) = c;
        else
            % 9: COMPRESSION: replace all simplex vertices
            X = delta.*X(1,:) + (1-delta).*X;            
        end
    end
    % delete from plot
    pause(1)
    delete(pS)
    tries = tries-1;
end

function plotObj = plotSimplex(X)
    plotObj = plot( [X(:,1);X(1,1)] , [X(:,2);X(1,2)] ,'r');
end