function res = reorderSinCos(osc)
%REORDERSINCOS reorders sine/cosine oscillation components
% input: oscillation components as vector
% output: oscillation components as array
% (number of coefficients must be even)
if mod(osc,2) == 1
    error('reorderSinCos: input error: oscillation component vector must contain even number of elements')
end
% S,C in that order 
res = [osc(1:2:end); osc(2:2:end)];
end