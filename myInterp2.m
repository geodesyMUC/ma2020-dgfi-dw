function out = myInterp2(x)
load('temp/output.mat');
tau1 = years(days(1:10:200));
tau2 = years(days(250:30:730));
out = interp2(tau1,tau2,resultGrid,x(:,1),x(:,2));
end