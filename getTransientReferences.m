function [tsT] = getTransientReferences(t_ts, types, low, upp, doStopTs)
%Creates Struct for Transients
% Input:
%   t_ts: vector with datetimes of events
%   types: cell with either 1 x nTransient or nEq x nTransient
%   low: array with lower constraints
%   upp: array with upper constraints
%   doStopTs: boolean, flag to stop transients from overlaying subsequent
%       eqs
% Output:
%   tsT: table with 4 columns {'time','type','lBound','uBound'}

if ~isempty(t_ts)
    if size(types, 1) == 1 && ...
            size(types, 2) == size(low, 2) && ...
            size(types, 2) == size(upp, 2)
        % type for all eq
        isGlobal = true;
    elseif length(t_ts) == size(types, 1) && ...
            size(types, 1) == size(low, 1) &&...
            size(types, 1) == size(upp, 1)
        % type for each eq
        isGlobal = false;
    elseif length(t_ts) ~= size(types, 1)
        error('getTransientReferences: number of transient types :  number of earthquakes mismatch')
    elseif size(types, 1) ~= size(upp, 1) || size(types, 1) ~= size(low, 1)
        error('getTransientReferences: number of constraints :  number of earthquakes mismatch')
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
            if doStopTs && i<length(t_ts)
                dt = t_ts(i+1) - t_ts(i);
                fprintf('getTransientReferences: input eq datetimes unit has to be YEARS: delta = %.3f\n', dt);
                fprintf('getTransientReferences: input transient lower limit unit has to be YEARS: %.3f\n', upp(i,j));
                fprintf('getTransientReferences: input transient upper limit unit has to be YEARS: %.3f \n', low(i,j));
                if low(i,j) > abs(dt)
                    % skip this transient
                    fprintf('getTransientReferences: transient %d for event %d removed\n', j, i);
                    continue
                elseif upp(i,j) > abs(dt)
                    % use dt as new constraint
                    fprintf('getTransientReferences: transient %d  for event %d adapted\n', j, i);
                    uppNew = abs(dt);
                end
            end
            % Append
            tsT = [tsT; { t_ts(i), types{eqidx,j}, low(eqidx,j), uppNew }];
        end
    end
end 
end

