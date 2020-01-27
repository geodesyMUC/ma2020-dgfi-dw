function [jumps] = getRelativeJumps(t, t0)
%getJumpsFromTable creates a jump vector using seconds since time series
% start (t0)
%   t: datetime vector with jumps for *current* station
%   t0: first observation in time series as a datetime

% subtract time first TS observation to get relative times [seconds]
jumps0 = etime(datevec(t), ...
    datevec(repmat(t0, size(t, 1), 1)));
% remove negative values -> jumps BEFORE Time Series starts
jumps0x = jumps0(jumps0 >= 0);
% sort jumps
jumps0x = sort(jumps0x);% ./(365.25 * 86400));

% resulting relative jumps vector %(can be converted to years if needed)
jumps = [jumps0x]; % jumps0x./(365.25 * 86400)

end

