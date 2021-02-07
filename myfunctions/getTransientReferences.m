function [tsT] = getTransientReferences(t_ts, types, low, upp, doRemoveTs)
%GETTRANSIENTREFERENCES Creates Struct for Transients
% Input:
%   t_ts: vector with datetimes of events
%   types: cell with either 1 x nTransient or nEq x nTransient
%   low: array with lower constraints
%   upp: array with upper constraints
%   doRemoveTS: transient will be removed if lower limit > time to next eq
%   doShiftLim: boolean, flag to prevent transients from overlaying subsequent
%       eqs by shifting the upper limit. 
%       Warning: For this to be applied, doRemoveTs also needs to be set to 
%       true. If not, it will be ignored.
%
% Output:
%   tsT: table with 4 columns {'time','type','lBound','uBound'}

tsT = cell2table(cell(0,4), 'VariableNames', {'time','type','lBound','uBound'});

% check out function input
if ~isempty(t_ts)
    if size(types, 1) == 1 && ...
            sum(~strcmp(types, '')) == size(low, 2) && ...
            sum(~strcmp(types, '')) == size(upp, 2)
        % input is type for ALL eq (default taus)
        fprintf('getTransientReferences: type of transient input recognized as global\n')
        isGlobal = true;
    elseif length(t_ts) == size(types, 1) && ...
            size(types, 1) == size(low, 1) &&...
            size(types, 1) == size(upp, 1)
        % input is type for EACH eq (custom taus)
        fprintf('getTransientReferences: type of transient input recognized as local\n')
        isGlobal = false;
    elseif any(~strcmp(types, '')) && size(types, 1) == 1 && ...
            size(types, 2) == size(low, 2) && ...
            size(low, 2) == size(upp, 2)
        % input transients < lower limit/upper limit. adjust limits 
        fprintf('getTransientReferences: number of transient types ~=  number of limits.')
        %low = low(~strcmp(types, ''));
        %upp = upp(~strcmp(types, ''));
        isGlobal = true;
    elseif ~any(~strcmp(types, ''))
         fprintf('getTransientReferences: no transients to be estimated\n')
         return
    elseif length(t_ts) ~= size(types, 1)
        error('getTransientReferences: input mismatch: number of transient types ~=  number of earthquakes') 
    elseif size(types, 1) ~= size(upp, 1) || size(types, 1) ~= size(low, 1)
        error('getTransientReferences: input mismatch: number of constraints ~=  number of earthquakes')
    end
else
    return
end
fprintf('getTransientReferences: global tau is "%s"\n', mat2str(isGlobal)); 
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
                fprintf('getTransientReferences: input transient lower limit unit must be YEARS. check: %.3f\n', upp(eqidx,j));
                fprintf('getTransientReferences: input transient upper limit unit must be YEARS. check: %.3f\n', low(eqidx,j));
                if doRemoveTs && low(eqidx,j)>abs(dt)
                    % skip this transient
                    fprintf('getTransientReferences: transient %d for event %d removed from model\n', j, i);
                    continue
                elseif doRemoveTs && upp(eqidx,j)>abs(dt)
                    warning('getTransientReferences:mismatch', ' transient %d  for event %d upper limit overwritten\n', j, i)
                    % use dt as new constraint ? set new value for uppNew
%                     uppNew = abs(dt);
%                     uppNew = 2*low(i,j); 
                end
            end
            % Append
            tsT = [tsT; { t_ts(i) , types{eqidx,j}, low(eqidx,j), uppNew }];
        end
    end
end 
end

