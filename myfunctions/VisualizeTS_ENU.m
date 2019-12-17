function VisualizeTS_ENU(data, STATION_NAME, jumps, jumpLogical)
% Visualize ENU TS
% Needs data (table with t, E, N, U)
% and station name string
titleString = sprintf('Station: %s (n of Measr. = %d)', STATION_NAME, size(data, 1));
for i = 1:3
    subplot(3, 1, i)
    plot(data{:, 'date'}, data{:, i + 2}, '.')
    hold on
    ylabel(sprintf('%s [mm]', data.Properties.VariableNames{2 + i}))
    xlim([min(data{:, 'date'}) max(data{:, 'date'})])
    grid on
    if i == 1, title(titleString), end
    if i == 3, xlabel('Time t'), end
    

    ax = gca;
    y1 = ax.YLim(1); % axis MIN
    y2 = ax.YLim(2); % axis MAX
    % Plot Red Line for all jumps
    if nargin == 3
        for j = 1:length(jumps)
            plot([jumps(j); jumps(j)], [y1; y2], 'r')
        end
    end
    
    hold off
end
end
