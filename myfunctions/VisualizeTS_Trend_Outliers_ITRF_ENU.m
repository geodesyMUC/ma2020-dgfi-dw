function VisualizeTS_Trend_Outliers_ITRF_ENU(t, ENU, ttrend, ENUtrend, titleString, outlierLogical, jumpTable, jumpITRF)
% Visualize ENU TS and Trend
% - date vector t
% - data (array with E, N, U)
% - trend (array with E, N, U)
% - string for plot titles, 1x3 str cell (E,N,U)
ENUstring = {'E', 'N', 'U'};

for i = 1:3
    subplot(3, 1, i)
    % Plot Measurements
    pPts = plot(t(~outlierLogical{i}), ENU(~outlierLogical{i}, i), ...
        '.', 'markersize', 4, ...
        'Color', [0, 0.4470, 0.7410]);
    hold on

    % Plot Trend
    pTrend = plot(ttrend, ENUtrend(:, i), 'r');
    ylabel(sprintf('%s [mm]', ENUstring{i}))
    xlim([min(t) max(t)])
    title(titleString{i});
    
    % Plot Outliers if there are any
    if nnz(outlierLogical{i}) > 0
        % Outlier
        pOutl = plot(t(outlierLogical{i}), ENU(outlierLogical{i}, i), ...
            '.', 'markersize', 8, ...
            'Color', [255, 153, 0]./255);
    end
    
    ax = gca;
    y1 = ax.YLim(1); % axis MIN
    y2 = ax.YLim(2); % axis MAX
    
    grid on
    hold on
    % if i == 3, xlabel('Time t'), end % only below 3rd plot
    
    % plot itrf jump vertical lines
    for j = 1:length(jumpITRF)
        pITRF = plot([jumpITRF(j); jumpITRF(j)], [y1; y2], '--', ...
            'color', [160, 160, 200]./255);
    end
    
    % Start Building legend entries
    lgdElements = [pPts, pTrend, pITRF];
    
    % Plot Jump Types
    pLogical = false(7,1); % Plot Logical for legend
    
    % Observations, 1:2 correspondends to obs and trend, 3:4 to itrf
    % datetime
    pLogical(1:3) = true; 
    
    % Plot Jumps from Jump Table
    for j = 1:size(jumpTable, 1)
        if jumpTable{j, 4} > 0
            % Earthquake
            pEq = plot([jumpTable{j, 2}; jumpTable{j, 2}], [y1; y2], ...
                'color', [102, 51, 0]./255);
            pLogical(5) = true;
        elseif jumpTable{j, 5} > 0
            % Antenna Change, ...
            pAnt = plot([jumpTable{j, 2}; jumpTable{j, 2}], [y1; y2], ...
                'color', [0.75, 0, 0.75]);
            pLogical(6) = true;
        else
            % Unknown cause
            pUnk = plot([jumpTable{j, 2}; jumpTable{j, 2}], [y1; y2], ...
                'color', [0.4660, 0.6740, 0.1880]);
            pLogical(7) = true;
        end
    end
    
    % Update legend Items Outliers
    if nnz(outlierLogical{i}) > 0
        % Outlier
        pLogical(4) = true;
        lgdElements = [lgdElements, pOutl];
    end
    
    % Update legend Items Jumps
    if pLogical(5) % eq
        lgdElements = [lgdElements, pEq];
    end
    if pLogical(6) % hw
        lgdElements = [lgdElements, pAnt];
    end
    if pLogical(7) % unk
        lgdElements = [lgdElements, pUnk];
    end
    
    % Item names default cell
    pNames = {'Observation', 'Trend', 'new ITRF', 'Outlier', 'Earthquake', 'HW Change', 'Unk. Cause'};
    % Create Legend
    legend(lgdElements, pNames(pLogical), 'Location', 'eastoutside')
    
    hold off
end
end
