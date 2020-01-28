function [] = writeOutputLog(fID, ts_name, xEst, polynDeg, P, HJumps, EQJump, rms, wrms)
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

% get count of unknowns for each type (poly, osc, jumps, transients)
nPolynTerms = polynDeg + 1; % 0, 1, 2, ... 
nOscParam = length(P) * 2; % cos & sin components (C, S) for every oscillation
nJumps = length(HJumps); % All Jumps - From DB and ITRF (if set to true)
nEQJumps = length(EQJump); % Only EQ Jumps -> n of transients

fprintf(fID, '\n### OUTPUT RESULTS FOR "%s" COORDINATE, ITERATIVELY REWEIGHTED LEAST SQUARES ###\n', ts_name);

strP = sprintf('p(%d): %.5f | ', [1:nPolynTerms; xEst(1:nPolynTerms)']); % Polynomial Coefficients
if ~isempty(P)
    strW = sprintf('w(%d): A = % .2fmm, C = % .2fmm,  S = % .2fmm\n', ...
        [1:length(P); ...
        sqrt(xEst(nPolynTerms + 1:2:nPolynTerms + nOscParam)'.^2 + ...
        xEst(nPolynTerms + 2:2:nPolynTerms + nOscParam)'.^2); % Amplitude A
        xEst(nPolynTerms + 1:2:nPolynTerms + nOscParam)'; ... % Oscillations C
        xEst(nPolynTerms + 2:2:nPolynTerms + nOscParam)']); % Oscillations S
else
    strW = 'None estimated';
end
if ~isempty(HJumps)
    strH = sprintf('J(%d): % .2fmm\n', [1:nJumps; ...
        xEst(nPolynTerms + nOscParam + 1:nPolynTerms + nOscParam + nJumps)']); % Jumps
else
    strH = 'None estimated';
end
if ~isempty(EQJump)
    strEQ = sprintf('A(%d): % .2fmm\n', [1:nEQJumps; ...
        xEst(...
        nPolynTerms + nOscParam + nJumps + 1:...
        nPolynTerms + nOscParam + nJumps + nEQJumps)']); % EQ Jumps
else
    strEQ = 'None estimated';
end
% Print all Parameters to log file
fprintf(fID, ['Estimated Parameters:\nPolynomial Coefficients:\n%s\n', ...
    'Oscillation Amplitudes:\n%s\nHeaviside Jumps:\n%s\nLogarithmic Transients:\n%s\n', ...
    'Root Mean Square Error RMS           = %.3fmm\n', ...
    'Weighted Root Mean Square Error WRMS = %.3fmm\n\n'], ...
    strP, strW, strH, strEQ, rms, wrms);
end

