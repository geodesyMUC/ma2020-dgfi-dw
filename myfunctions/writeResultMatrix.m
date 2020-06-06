function [resultM] = writeResultMatrix(t_trend,trends, doITRFjump, KK, p, outl_factor, rms_enu, wrms_enu)
%writeResultMatrix writes results of a 3fold IRLS Trend Estimation (ENU/XYZ) to a csv file
if size(t_trend, 1) < size(t_trend, 2) 
    t_trend = t_trend'; % transpose 
end
% Set up matrix with results and LSE parameters
resultM = [posixtime(t_trend),  trends];
resultM = [[200, wrms_enu(1), wrms_enu(2), wrms_enu(3)]; ...
    resultM]; % 200: WRMS (third line)
resultM = [[100, rms_enu(1), rms_enu(2), rms_enu(3)]; ...
    resultM]; % 100: RMS (second line)
resultM = [[NaN, KK, p, outl_factor]; resultM]; % LSE parameters (first line)
% NaN used to be doITRF in previous version, needs fixing since doITRF is
% vector now
end

