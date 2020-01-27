function [itrf_jumps] = getRelativeITRFJumps(t0, itrf_txt_path)
%getJumpsFromTable creates a jump vector using seconds since time series
% start (t0)
%   t0: first observation in time series as a datetime
%   itrf_txt_path: path to simple text file containing itrf change dates as
%   string (yyyy-MM-dd UTC)
%   This function needs the function "getRelativeJumps" to calculate relative datetimes.

% ITRF Jumps 
ditrf = readITRFChanges(itrf_txt_path);
% Apply getRelativeJumps function to get relative time to t0 of TS in
% [seconds]
itrf_jumps = getRelativeJumps(ditrf, t0);

end

