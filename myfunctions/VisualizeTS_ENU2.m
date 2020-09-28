function VisualizeTS_ENU2(data, station, jumps, varargin)
% Visualize ENU TS
% data, STATION_NAME, jumps, jumpLogical
% Needs data (table with t, E, N, U)
% and station name string

% Defaults
defaultJumpsTypes = zeros(size(jumps, 1), 2);

% Input Parser Object
p = inputParser;
% Validate Arguments
validJumpsTypesArray = @(x) size(x, 1) == size(jumps,1) && size(x, 2) == 2 && ...
    isnumeric(x);

addRequired(p, 'data', @(x) size(x, 2)==5);
addRequired(p, 'station', @(x) ischar(x));
addRequired(p, 'jumps', @(x) isdatetime(x));

addParameter(p, 'jumpTypes', defaultJumpsTypes, validJumpsTypesArray);

% Parse Input Parameters
parse(p, data, station, jumps, varargin{:});

titleString = sprintf('TS Observations, Station: %s (n of Measr. = %d)', p.Results.station, size(p.Results.data, 1));

for i = 1:3
    subplot(3, 1, i)
    pPts = plot(p.Results.data{:, 'date'}, p.Results.data{:, i + 2}, '.', ...
        'Color', [0, 0.4470, 0.7410]); % blue
    hold on
    ylabel(sprintf('%s [mm]', p.Results.data.Properties.VariableNames{2 + i}))
    xlim([min(p.Results.data{:, 'date'}) max(p.Results.data{:, 'date'})]) % red
    grid on
    if i == 1, title(titleString), end
    if i == 3, xlabel('Time t'), end
    
    ax = gca;
    y1 = ax.YLim(1); % axis MIN
    y2 = ax.YLim(2); % axis MAX
    
    % Plot Jump Types
    pLogical = false(4,1); % Plot Logical for legend
    pLogical(1) = true; % Observations
    for j = 1:size(p.Results.jumpTypes, 1)
        if p.Results.jumpTypes(j, 1) > 0
            % Earthquake
            pEq = plot([p.Results.jumps(j); p.Results.jumps(j)], [y1; y2], ...
                'color', [1, 0, 0]);
            pLogical(2) = true;
        elseif p.Results.jumpTypes(j, 2) > 0
            % Antenna Change, ...
            pAnt = plot([p.Results.jumps(j); p.Results.jumps(j)], [y1; y2], ...
                'color', [0.75, 0, 0.75]);
            pLogical(3) = true;
        elseif ~isnat(p.Results.jumps(j))
            % Unknown cause
            pUnk = plot([p.Results.jumps(j); p.Results.jumps(j)], [y1; y2], ...
                'color', [0.4660, 0.6740, 0.1880]);
            pLogical(4) = true;
        end
    end
    hold off
    
    % Build legend entries
    lgdElements = [pPts];
    
    if pLogical(2) == 1 % eq
        lgdElements = [lgdElements, pEq];
    end
    if pLogical(3) == 1 % hw
        lgdElements = [lgdElements, pAnt];
    end
    if pLogical(4) == 1 % unk
        lgdElements = [lgdElements, pUnk];
    end   
    pNames = {'Observation', 'Earthquake', 'HW Change', 'Unk. Cause'};
    legend(lgdElements, pNames(pLogical), 'Location', 'eastoutside')
end
end
