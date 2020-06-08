function [] = writeInputLog(fID, stationname, coordinateName, t, polynDeg, P, HJumps, t_trans, taus, transType, KK, p, outl_factor)
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
t_end = t(end);

% header
fprintf(fID, '### INPUT: Station %s, Coordinate "%s" ###\n', stationname, coordinateName);
fprintf(fID, 'Time Series Start t_0  = %s\n', datestr(t0, 'yyyy-mm-dd HH:MM'));
fprintf(fID, 'Time Series End t_end  = %s\n\n', datestr(t_end, 'yyyy-mm-dd HH:MM'));
% polynomial model
if polynDeg<0; polynStr='none estimated'; else; polynStr=sprintf('%d', polynDeg); end
fprintf(fID, 'Degree of Polynomial:\n%s\n\n', polynStr) ;
% oscillation model
if isempty(P); oscStr=['none estimated', '\n']; else; ...
        oscStr=sprintf('%.1fy | %.2fd\n', [P; P.*365.25]); end
fprintf(fID, 'Oscillations:\n%s\n', oscStr);
% heaviside jump model
fprintf(fID, 'Heaviside Jumps [Relative Time - Datetime]:');
if isempty(HJumps); fprintf(fID, '\nnone estimated');end
fprintf(fID, '\n');
for i = 1:length(HJumps)
    % Print Jump Datetime Information
    fprintf(fID, 't0 + %.2fs | %s\n', HJumps(i), ...
        datestr(t0 + seconds(HJumps(i)), 'yyyy-mm-dd HH:MM'));
end
% transient model - datetimes
fprintf(fID, '\nTransients [Relative Time - Datetime - Type of Transient]:');
if isempty(t_trans); fprintf(fID,'\nnone estimated');end
for i = 1:length(t_trans)
    transPrint = cellstr(transType);
    % Print Jump Datetime Information
    fprintf(fID, '\nt0 + %.2fs | %s | %s ', t_trans(i), ...
        datestr(t0 + seconds(t_trans(i)), 'yyyy-mm-dd HH:MM'), ...
        sprintf('%s ', transPrint{:})...
        );
end

% transient model - tau values
fprintf(fID, '\n\nTransient 1 Range [Lower Bound | Upper Bound]:\n'); % 1st Transient
if isempty(taus)
    fprintf(fID,'not estimated\n');
else
    fprintf(fID, '%dd | %dd\n', days(years(min(taus(:,1)))), days(years(max(taus(:,1)))));
end
fprintf(fID, '\nTransient 2 Range [Lower Bound | Upper Bound]:\n');   % 2nd Transient
if size(taus,2)<2
    fprintf(fID,'not estimated\n');
else
    fprintf(fID, '%dd | %dd\n', days(years(min(taus(:,2)))), days(years(max(taus(:,2)))));
end

% LSE/IRLSE
fprintf(fID, ['\nIRLSE/LSE Parameters:\nKK = %d (n of Iterations in IRLSE)\np = %.1f (L_p Norm used for IRLS)\n', ...
    'outl_factor = %d (median of error + standard deviation * factor < outlier)\n\n'], KK, p, outl_factor);
end

