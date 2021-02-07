function y = TimeFunction(t, p, cs, w, j_t, b, ts_t, a, tsTau, tsType, doTsOverlay)
% Creates time series from estimated parameters
% t: timestamps [years]
% polyC: polynomial coeffiecients
% cs: cos/sin periodic coefficients
% w: vector containing periods (rad)
% shift
% b: jump terms
% eq jump times [years]
% amplitude of transient
% relaxation time T [years]
% type of tau ("log"|"exp")

if size(t, 1) > size(t, 2) % fix row/col 
    t = t';
end
p = fixRowToCol(p);
cs = fixRowToCol(cs);
b = fixRowToCol(b);
a = fixRowToCol(a);

nPoly = length(p); % number of polynomial coeffiecients
nW = length(w); % number of periodic coefficients
nJump = length(j_t); % number of shifts
nTs = length(ts_t); % number of transients

if length(tsTau) ~= nTs || nTs ~= length(tsType)
    error('Time Function:transient model mismatch: vector length of tau,tau datetime and type of tau')
end

% test
A_ = createCoeffMat(t, nPoly-1, w, j_t, ts_t, tsTau, tsType, doTsOverlay);
x = [p, cs, b, a];
y = A_*x';

end

function in = fixRowToCol(in)
if size(in, 1)>size(in, 2)
    in = in';
end
end