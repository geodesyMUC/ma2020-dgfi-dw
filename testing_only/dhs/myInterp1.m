function out = myInterp1(x)
load('src/output.mat');
tau1 = years(days(1:10:200));
tau2 = 15;
resultRow = resultGrid(tau2,:);
out = interp1(tau1,resultRow,x);
end