function [jumps_eq] = getRelativeJumps_eq(t, t0, eq_logical)
%getJumpsFromTable_eq creates an earthquake jump vector using seconds since time series
% start (t0)
%   t: datetime vector with jumps for *current* station
%   t0: datetime of first observation in time series
%   eq_logical: logical where 1:= earthquake and ~1:=no earthquake
%
%   This function needs the function "getRelativeJumps" to calculate relative datetimes.

if nargin == 2
    % all jumps are eq
    jumps_eq = getRelativeJumps(t, t0);
elseif nargin == 3
    % logical denotes eq jumps, others will be ignored
    jumps_eq = getRelativeJumps(t(eq_logical), t0);
end

end
