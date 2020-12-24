function [] = writeOutputLog(fID, station, cnName, polyn, osc, jump, t_ts, tsAmp, tsTau, outl, errors)
%writeOutputLog writes output parameters to a formatted log file
%   fID: log file identifier
%   ts_name: Name of Time Series
%   xEst: estimated parameters
%   polynDeg: integer power of polynome denoting station velocity
%   P: vector containing periods in [YEARS]
%   HJumps: vector containing jump times in [SECONDS] relative to t0 (if not specified -> empty)
%   EQJump: vector containing jump times for earthquakes in [SECONDS] relative to t0 (if not specified -> empty)
%   erros: key-value cell
fprintf(fID, '## RESULTS: Station %s, Coordinate "%s" -----------------------------------\n', station, cnName);

if ~isempty(polyn)
    strP = sprintf('%.2f\n', polyn);
else
    strP = sprintf('%s\n', 'none estimated');
end
fprintf(fID, '# Polynomial Coefficients\n%s\n', strP);

if ~isempty(osc)
    osc = reorderSinCos( osc ); % S*cos(wt), C*sin(wt)
    strO = '';
    for i = 1:size(osc,2)
        strONew = sprintf('%.2f; %.2f; %.2f; %.2f\n', osc(:,i), ...
            sqrt(osc(1,i)^2+osc(2,i)^2), ...
            atan( osc(1,i)/osc(2,i) )  );
        strO = [strO strONew];
    end
else
    strO = sprintf('%s\n', 'none estimated');
end
fprintf(fID, '# Oscillation Amplitude Components (S; C; A; phase shift)\n%s\n', strO);

if ~isempty(jump)
    strJ = sprintf('%.2f\n', jump);
else
    strJ = sprintf('%s\n', 'none estimated', '\n');
end
fprintf(fID, '# Heaviside Jump Coefficients\n%s\n', strJ);

% reshape tau/ts amp
nTs = length(t_ts); % number of transients

% tau relaxation times
if ~isempty(tsTau)
    strR = '';
    for i = 1:nTs
        strRnew = sprintf('%.2f\n', days(years( tsTau(i) )));
        strR = [strR, strRnew];
    end
else
    strR = sprintf('%s\n', 'none estimated');
end
fprintf(fID, '# Transient Relaxation Times Tau (days)\n%s\n', strR);

% tau amplitudes
if ~isempty(tsAmp)
    strT = '';
    for i = 1:nTs
        strTNew = sprintf('%.2f\n', tsAmp(i) );
        strT = [strT, strTNew];
    end
else
    strT = sprintf('%s\n', 'none estimated');
end
fprintf(fID, '# Transient Amplitudes\n%s\n', strT);

% outliers
strOutl = sprintf('%d\n', sum(outl) );
fprintf(fID, '# Outliers (number of measurements removed)\n%s\n', strOutl);

% errors
fprintf(fID, '# Error metrics\n');
for i = 1:size(errors,1)
    fprintf(fID, '%s: %.3f\n', pad(errors{i,1}, 8),errors{i,2});
end
fprintf(fID, '\n');
end

