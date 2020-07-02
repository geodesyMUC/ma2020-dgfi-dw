function [tsT] = getTransientReferences(t_ts, types, low, upp, doShiftLim)
%Creates Struct for Transients
% Input:
%   t_ts: vector with datetimes of events
%   types: cell with either 1 x nTransient or nEq x nTransient
%   low: array with lower constraints
%   upp: array with upper constraints
%   doShiftLim: boolean, flag to prevent transients from overlaying subsequent
%       eqs by shifting the upper limit
% Output:
%   tsT: table with 4 columns {'time','type','lBound','uBound'}

% check out function input
if ~isempty(t_ts)
    if size(types, 1) == 1 && ...
            size(types, 2) == size(low, 2) && ...
            size(types, 2) == size(upp, 2)
        % input is type for ALL eq (default taus)
        fprintf('getTransientReferences: type of transient input recognized as global\n')
        isGlobal = true;
    elseif length(t_ts) == size(types, 1) && ...
            size(types, 1) == size(low, 1) &&...
            size(types, 1) == size(upp, 1)
        % input is type for EACH eq (custom taus)
        fprintf('getTransientReferences: type of transient input recognized as local\n')
        isGlobal = false;
    elseif length(t_ts) ~= size(types, 1)
        error('getTransientReferences: input mismatch: number of transient types ~=  number of earthquakes')
    elseif size(types, 1) ~= size(upp, 1) || size(types, 1) ~= size(low, 1)
        error('getTransientReferences: input mismatch: number of constraints ~=  number of earthquakes')
    end
end

tsT = cell2table(cell(0,4), 'VariableNames', {'time','type','lBound','uBound'});
for i = 1:length(t_ts)
    if isGlobal
        eqidx = 1;
    else
        eqidx = i;
    end
    for j = 1:size(types,2)
        if ~isempty( types{eqidx,j} )
            uppNew = upp(eqidx,j);
            % Adjust constraints if flag is set
            if i<length(t_ts)
                dt = t_ts(i+1) - t_ts(i);
                fprintf('getTransientReferences: input eq datetimes unit must be YEARS. check: delta t = %.3f\n', dt);
                fprintf('getTransientReferences: input transient lower limit unit must be YEARS. check: %.3f\n', upp(i,j));
                fprintf('getTransientReferences: input transient upper limit unit must be YEARS. check: %.3f\n', low(i,j));
                if low(i,j)>abs(dt)
                    % skip this transient
                    fprintf('getTransientReferences: transient %d for event %d removed from model\n', j, i);
                    continue
                elseif doShiftLim && upp(i,j)>abs(dt)
                    % use dt as new constraint
                    fprintf('getTransientReferences: transient %d  for event %d upper limit overwritten\n', j, i);
                    uppNew = abs(dt);
%                     uppNew = low(i,j)+1e-3;
                end
            end
            % Append
            tsT = [tsT; { t_ts(i) , types{eqidx,j}, low(eqidx,j), uppNew }];
        end
    end
end 
end

