function [] = writeInputLog(fID, stationname, coordinateName, t, polynDeg, P, HJumps, EQJump, KK, p, outl_factor)
%writeInputLog writes input parameters to a formatted log file
%   fID: log file identifier
%   stationname: name of the station
%   coordinateName: name of coordinate
%   t: time series observation datetime vector
%   polynDeg: integer power of polynome denoting station velocity
%   P: vector containing periods in [YEARS]
%   HJumps: vector containing jump times in [SECONDS] relative to t0 (if not specified -> empty)
%   EQJump: vector containing jump times for earthquakes in [SECONDS] relative to t0 (if not specified -> empty)
%   KK: number of iterations for IRLS
%   p: L_p Norm for IRLS. p=2 equals to the euclidean norm (=L2). 
%   If p=2, no reweighting will be applied, independent of the number of iterations KK
%   outl_factor: median(error)|mean(error) + standard deviation * factor -> outlier

t0   = t(1);
tend = t(end);

fprintf(fID, 'Time Series Analysis: Trend for Station "%s", Coordinate "%s" -----------\n', stationname, coordinateName);
fprintf(fID, 'Time Series Start t_0  = %s\n', datestr(t0, 'yyyy-mm-dd HH:MM'));
fprintf(fID, 'Time Series End t_end  = %s\n\n', datestr(tend, 'yyyy-mm-dd HH:MM'));

fprintf(fID, '\n### INPUT PARAMETERS ###\n');

fprintf(fID, 'Polynomial Trend Model (SLTM):\nDegree = %d\n\n', polynDeg) ;
fprintf(fID, 'Oscillations:\n%s\n\n', ...
    sprintf('w(%d) = %.1f y | %.2f d\n', [(1:length(P)); P; P.*365.25]));
% All jumps
fprintf(fID, 'Jump table (Heaviside):\n');
for i = 1:length(HJumps)
    % Print Jump Datetime Information
    fprintf(fID, 'J(%d) = t0 + %.2fs | %s\n', i, HJumps(i), ...
        datestr(t0 + seconds(HJumps(i)), 'yyyy-mm-dd HH:MM'));
end

% EQ
fprintf(fID, '\nEarthquake Jump Table for Heaviside (Offset) & Logarithmic Transient (Postseismic Movement) estimation:\n');
for i = 1:length(EQJump)
    % Print Jump Datetime Information
    fprintf(fID, 'J(%d) = t0 + %.2fs | %s\n', i, EQJump(i), ...
        datestr(t0 + seconds(EQJump(i)), 'yyyy-mm-dd HH:MM'));
end
% LSE/IRLSE
fprintf(fID, ['\nIRLSE/LSE Parameters:\nKK = %d (n of Iterations in IRLSE)\np = %.1f (L_p Norm used for IRLS)\n', ...
    'outl_factor = %d (median of error + standard deviation * factor < outlier)\n'], KK, p, outl_factor);
end

