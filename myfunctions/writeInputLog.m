function [] = writeInputLog(fID, stationname, coordinateName, t, polynDeg, P, HJumps, t_trans, lBoundTau, uBoundTau, transType, KK, p, outl_factor)
%writeInputLog writes input parameters to a formatted log file
%   fID: log file identifier
%   stationname: name of the station
%   coordinateName: name of coordinate
%   t: time series observation datetime vector
%   polynDeg: integer power of polynome denoting station velocity
%   P: vector containing periods in [YEARS]
%   HJumps: vector containing jump times in [YEARS] relative to t0 (if not specified -> empty)
%   EQJump: vector containing jump times for earthquakes in [YEARS] relative to t0 (if not specified -> empty)
%   KK: number of iterations for IRLS
%   p: L_p Norm for IRLS. p=2 equals to the euclidean norm (=L2). 
%   If p=2, no reweighting will be applied, independent of the number of iterations KK
%   outl_factor: median(error)|mean(error) + standard deviation * factor -> outlier

t0   = t(1);
t_end = t(end);

% header
fprintf(fID, '## INPUT MODEL: Station %s, Coordinate "%s" ================================\n', stationname, coordinateName);
% fprintf(fID, 'Time Series Start t_0: %s\n', datestr(t0, 'yyyy-mm-dd HH:MM'));
% fprintf(fID, 'Time Series End t_end: %s\n\n', datestr(t_end, 'yyyy-mm-dd HH:MM'));

% polynomial model
if polynDeg<0; polynStr='none estimated'; else; polynStr=sprintf('%d', polynDeg); end
fprintf(fID, '# Degree of Polynomial\n%s\n\n', polynStr) ;

% oscillation model
if isempty(P)
    strW = sprintf('%s\n', 'none estimated');
else
    strW = sprintf('%.1f; %.2f\n', [P.*365.25; P]);
end
fprintf(fID, '# Oscillations [days; years]\n%s\n', strW);

% heaviside jump model
if isempty(HJumps)
    strJ = sprintf('%s\n', 'none estimated');
else
    strJ = '';
    for i = 1:length(HJumps)
        % Print Jump Datetime Information
        strJNew = sprintf('%s; %.4f\n', ...
            datestr(t0 + years(HJumps(i)), 'yyyy-mm-dd HH:MM'), ...
            HJumps(i));
        strJ = [strJ, strJNew];
    end
end
fprintf(fID, '# Heaviside Jumps (Rel. to t0 [years];  Datetime])\n%s\n', strJ);

% transient model - datetimes
if isempty(t_trans)
    strTs = sprintf('%s\n', 'none estimated');
else
    strTs = '';
    transPrint = cellstr(transType);
    for i = 1:length(t_trans)
        % Print Jump Datetime Information
        strTsNew = sprintf('%s; %.4f; %s \n', ...
            datestr(t0 + years(t_trans(i)), 'yyyy-mm-dd HH:MM'), t_trans(i), ...
            sprintf('%s ', transPrint{i}) );
        strTs = [strTs, strTsNew];        
    end
end
fprintf(fID, '# Transients (rel. to t0 [years]; Datetime; Type of Transient)\n%s\n', strTs);

% transient model - tau value range
if isempty(t_trans)
    strTsB = sprintf('%s\n', 'none estimated');
else
    strTsB = '';
    for i = 1:length(t_trans)
        strTsBNew = sprintf('%.1f; %2.4f; %.1f; %2.5f\n', ...
            days(years( lBoundTau(i) )), lBoundTau(i),  ...
            days(years( uBoundTau(i) )), uBoundTau(i));
        strTsB = [strTsB, strTsBNew];
    end
end
fprintf(fID, '# Transient Range (Lower Limit [days;years] ; Upper Limit [days;years])\n%s\n', strTsB); % 1st Transient

% LSE/IRLSE
fprintf(fID, ['# IRLSE/LSE Parameters\nKK: %d\nL_p Norm: %.1f\n', ...
    'outl_factor: %d\n\n'], KK, p, outl_factor);
end

