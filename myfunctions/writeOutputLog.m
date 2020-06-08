function [] = writeOutputLog(fID, station, cnName, polyn, osc, jump, ts, tau, rms, wrms)
%writeOutputLog writes output parameters to a formatted log file
%   fID: log file identifier
%   ts_name: Name of Time Series
%   xEst: estimated parameters
%   polynDeg: integer power of polynome denoting station velocity
%   P: vector containing periods in [YEARS]
%   HJumps: vector containing jump times in [SECONDS] relative to t0 (if not specified -> empty)
%   EQJump: vector containing jump times for earthquakes in [SECONDS] relative to t0 (if not specified -> empty)
%   rms: estimated root mean square
%   wrms: estimated weighted root mean square
fprintf(fID, '### RESULTS: Station %s, Coordinate "%s" ###\n', station, cnName);

if ~isempty(polyn)
    strP = sprintf('%.3f; ', polyn);
else
    strP = 'none estimated';
end
fprintf(fID, 'Polynomial Coefficients:\n%s\n\n', strP);

if ~isempty(osc)
    strO = '';
    for i = 1:size(osc,2)
        strONew = sprintf('%.3f; %.3f\n', osc(:,i));
        strO = [strO strONew];
    end
else
    strO = 'none estimated';
end
fprintf(fID, 'Oscillation Amplitudes [cosine sine]:\n%s\n', strO);

if ~isempty(jump)
    strJ = sprintf('%.3f;\n', jump);
else
    strJ = ['none estimated', '\n'];
end
fprintf(fID, 'Heaviside Jump Coefficients:\n%s\n', strJ);

if ~isempty(tau)
    strR = sprintf('%.1f;\n', days(years(tau)));
else
    strR = 'none estimated';
end
fprintf(fID, 'Transient Relaxation Times Tau:\n%s\n', strR);

nEq = size(ts,2)/size(tau,2); % number of eq
if ~isempty(ts)
    strT = '';
    ts = reshape(ts, [nEq size(tau,2)]);
    for i = 1:nEq
        strTNew = sprintf('%s\n', sprintf('%.3f; ', ts(i,:) ));
        strT = [strT, strTNew];
    end
else
    strT = 'none estimated';
end
fprintf(fID, 'Transient Amplitudes [tau1 tau2](per Earthquake):\n%s\n', strT);

fprintf(fID, 'RMS:\n%.3f;\nWRMS:\n%.3f;\n\n', rms, wrms);
end

