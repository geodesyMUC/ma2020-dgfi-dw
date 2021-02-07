function [] = writeHeaderLog(fID, stationname, t0, tn, n, estOpt, tf, doitrf, doTsOverlay, doRemoveTs)

fprintf(fID,'## Header ================================================================================\n');

strSt = sprintf('%s\n', stationname);
fprintf(fID,'# Station Name\n%s\n', strSt);

strTs = sprintf('%s; %s\n', datestr(t0, 'yyyy-mm-dd HH:MM:SS'), datestr(tn, 'yyyy-mm-dd HH:MM:SS') );
fprintf(fID,'# Time Series (Start;End)\n%s\n', strTs);

strObs = sprintf('%d\n', n);
fprintf(fID,'# Number of Measurements\n%s\n', strObs);

strEst = sprintf('%d\n', estOpt);
fprintf(fID,['# LS Estimation Method (1 := E,N,U separately;'...
    ' 2 := E,N combined, U sep.; 3 := E,N,U comb.)\n%s\n'], strEst);

strTar = sprintf('%s\n', tf);
fprintf(fID,'# Target Function for Optimization\n%s\n', strTar);

strItrf = sprintf('%s; %s; %s\n', ...
    mat2str( doitrf(1) ), mat2str( doitrf(2) ), mat2str( doitrf(3) ) );
fprintf(fID,'# new ITRF jump\n%s\n', strItrf);

strTsopt = sprintf('doTsOverlay: %s\ndoRemoveTs: %s\n', mat2str(doTsOverlay), mat2str(doRemoveTs));
fprintf(fID,'# Transient Settings\n%s\n', strTsopt);

fprintf(fID,'# Weighting Settings (global)\n%s\n', 'No weighting applied');

end

 