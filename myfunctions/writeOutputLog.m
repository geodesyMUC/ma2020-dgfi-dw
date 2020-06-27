function [] = writeOutputLog(fID, station, cnName, polyn, osc, jump, t_ts, tauType, tsAmp, tsTau, rms, wrms)
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
    strP = sprintf('%.3f;\n', polyn);
else
    strP = sprintf('%s\n', 'none estimated');
end
fprintf(fID, 'Polynomial Coefficients:\n%s\n', strP);

if ~isempty(osc)
    strO = '';
    for i = 1:size(osc,2)
        strONew = sprintf('%.3f; %.3f\n', osc(:,i));
        strO = [strO strONew];
    end
else
    strO = sprintf('%s\n', 'none estimated');
end
fprintf(fID, 'Oscillation Amplitudes [cosine sine]:\n%s\n', strO);

if ~isempty(jump)
    strJ = sprintf('%.3f;\n', jump);
else
    strJ = sprintf('%s\n', 'none estimated', '\n');
end
fprintf(fID, 'Heaviside Jump Coefficients:\n%s\n', strJ);

% reshape tau/ts amp
nEq = length(t_ts); % number of eq
nAmp = nnz(tsAmp);
nTau = nnz(tsTau);
nTauType = size(tauType,1);
if nTau == nTauType
    strTau = 'global';
    tsAmp = reshape(tsAmp, [nEq size(tsTau,2)]);
    if nTauType > 1
        tsAmp = tsAmp';
    end
else
    strTau = 'local';
    tsTau = reshape(tsTau, [nEq nTauType]);
    tsAmp = reshape(tsAmp, [nEq nTauType]);
    if nTauType > 1
        tsAmp = tsAmp';
        tsTau = tsTau';
    end
end

% tau relaxation times
if ~isempty(tsTau) && strcmp(strTau, 'local')
    strR = '';
    for i = 1:nEq
        strRnew = sprintf('%s\n', sprintf('%.1f; ', days(years( tsTau(i,:) ))));
        strR = [strR, strRnew];
    end
elseif ~isempty(tsTau) && strcmp(strTau, 'global')
       strR = sprintf('%s\n', sprintf('%.1f; ', days(years( tsTau ))));
else
    strR = sprintf('%s\n', 'none estimated');
end
fprintf(fID, 'Transient Relaxation Times Tau (days), estimation %s:\n%s\n', strTau, strR);

% tau amplitudes
if ~isempty(tsAmp)
    strT = '';
    for i = 1:nEq
        strTNew = sprintf('%s\n', sprintf('%.3f; ', tsAmp(i,:) ));
        strT = [strT, strTNew];
    end
else
    strT = sprintf('%s\n', 'none estimated');
end

fprintf(fID, 'Transient Amplitudes [tau1 tau2](per Earthquake):\n%s\n', strT);

fprintf(fID, 'RMS:\n%.3f;\nWRMS:\n%.3f;\n\n', rms, wrms);
end

