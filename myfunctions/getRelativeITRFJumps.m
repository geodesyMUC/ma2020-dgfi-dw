function [itrf_jumps] = getRelativeITRFJumps(t0, itrf_txt_path)
%getJumpsFromTable creates a jump vector using seconds since time series
% start (t0)
%   t0: first observation in time series as a datetime
%

% ITRF Jumps 
ditrf = readITRFChanges(itrf_txt_path);

% keep format "yyyy-MM-dd", time will be 0:00 am.

jumps0itrf = etime(datevec(ditrf), ...
    datevec(repmat(t0, size(ditrf, 1), 1)));
% remove negative values -> itrf jumps BEFORE Time Series starts
jumps0itrf = jumps0itrf(jumps0itrf >= 0);
% sort asc
jumps0itrf = sort(jumps0itrf);% ./(365.25 * 86400));
% result itrf jumps vector %(can be converted to years if needed)
itrf_jumps = [jumps0itrf]; % jumps0x./(365.25 * 86400)
end

