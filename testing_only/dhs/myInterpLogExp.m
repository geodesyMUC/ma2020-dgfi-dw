function out = myInterpLogExp(x)
load('src/output-logexp2.mat');
tau1 = years(days(1:2:50));
tau2 = years(days(51:5:365*2));
out = interp2(tau1,tau2,resultGrid,x(:,1),x(:,2));
end