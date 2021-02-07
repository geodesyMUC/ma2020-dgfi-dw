function [jumps] = getRelativeJumps(tJ, t0, tL)
%getJumpsFromTable creates a jump vector [seconds] since time series
% start (t0)
%   t: datetime vector with jumps for *current* station
%   t0: first observation in time series as a datetime
%   tL: last observation in time series as a datetime
% subtract time first TS observation to get relative times [seconds]
jumps0 = etime(datevec(tJ), ... % before = negative
    datevec(repmat(t0, size(tJ, 1), 1)));

jumps1 = etime(datevec( repmat(tL, size(tJ,1), 1) ) , ...
    datevec(tJ)); % after = negative

% remove negative values -> jumps BEFORE Time Series starts
% remove negative values -> jumps AFTER Time Series ends
jumps0x = jumps0(jumps0 >= -604800 & jumps1 > 0);

% sort jumps
jumps0x = sort(jumps0x);

% resulting relative jumps vector %(can be converted to years if needed)
jumps = [jumps0x]; 
jumps = seconds(jumps);
end