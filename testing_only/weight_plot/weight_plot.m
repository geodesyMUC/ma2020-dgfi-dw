% just a small fct to make the weight plot
% control weight decay after eq.
close all
clear variables

load('tx.mat')

coordinateName = ["U"];
stationName = "21729S007A04";

ts = [0 0]; % transient dummys
wFactor = [-0.5 0 1];  % [-1;1], -1 := strong decay ; 0 := linear decay ; 1 := no decay
twEq  = 0;      % time at which weighting takes eq weight "wEq" in [years], rel. time to ts
twNo = 2;       % time at which weighting takes default weight "wNo" in [years], rel. time to ts
wEq = [1 1 1];        % weight at time t_ts (=eq)
wNo = [0.5 0.2 0];        % default weight

w = zeros(length(tx),1);

figure % debug
hold all
for k = 1:length(wFactor)
    w(:) = wNo(k);
    tTs = unique( ts );
    
    
    %         subplot(3,1,i) % debug
    for j = 1:length( tTs ) % Loop Transients
        twSta = [twEq; wEq(k)];
        twMid = [twEq + (twNo-twEq)*0.5 ; wNo(k) + (wEq(k)-wNo(k))*0.5];
        
        if ( wEq(k)-wNo(k) ) < ( twNo-twEq )
            scale =  norm( twSta-twMid ) * ( wEq(k)-wNo(k) )/( twNo-twEq );
        else
            scale =  norm( twSta-twMid ) * ( twNo-twEq )/( wEq(k)-wNo(k) );
        end
        
        
        if wFactor(k) >= 1
            warning('weighting: wFactor adjusted to be < 1')
            wFactor(k) = 1-1e-5;
        elseif wFactor(k) <= -1
            warning('weighting: wFactor adjusted to be > -1')
            wFactor(k) = -1+1e-5;
        end
        % v1:
        %             V = twSta - twMid;
        %             nuV = [ V(2); -V(1) ] ./ norm( V ) * scale * wFactor(k);
        %             twPivot = twMid + nuV;
        % v2:
        V = [twNo; wEq(k)] - twMid;
        newV = V * wFactor(k);
        twPivot = twMid + newV;
        % ---
        dt = tx - tTs(j);
        wLine = [...
            twEq, twPivot(1), twNo; ...
            wEq(k), twPivot(2), wNo(k)];
        wx = interp1(wLine(1,:), wLine(2,:), dt, 'pchip', NaN);
        wxIdx = ~isnan(wx);
        w( wxIdx) = wx( wxIdx );
        
        
        hold on % debug
        %plot([twMid(1)+tTs(j) twPivot(1)+tTs(j)],[twMid(2) twPivot(2)], 'r'); % debug
        pc=plot(twSta(1)+tTs(j), twSta(2), 'gx'); % debug
%         pc=plot(twMid(1)+tTs(j), twMid(2), 'gx'); % debug
        pc=plot(twPivot(1),twPivot(2),'gx');
        pc=plot(twNo, wNo, 'gx'); % debug
    end
    p(k) = plot(tx,w(:), '.'); % debug
    
end

ylabel('weight') % debug
xlabel('t rel. to t_{0} [y]') % debug
title(sprintf('Weightings for "%s" (%s)', coordinateName, stationName)) % debug
legend([pc(1); p(:)],{'weight controls','weighting1','weighting2','weighting3'})