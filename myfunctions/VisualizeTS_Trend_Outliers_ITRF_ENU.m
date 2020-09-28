function VisualizeTS_Trend_Outliers_ITRF_ENU(t_obs, obs, obs_outlier_logical, obs_names, ...
    t_trend, trends, title_strings, t_jump, jump_categories_logical, jump_categories_names, jumpITRF)
%VisualizeTS_Trend_Outliers_ITRF_ENU: Visualizes Time Series Observations, calculated trends and specified jumps.
%   INPUT
%   (n x 1) t_obs: datetime vector of raw observations
%   (n x m) obs: raw observations. m denotes the number of different time
%           series, i.e. E, N and U.
%   (n x m) obs_outlier_logical: logical for time series where 1:=outlier
%           and 0:=no outlier.
%   {1 x m} obs_names: name strings of time series observations, i.e. East,
%           North, Up. This is used in plot legend.
%   (1 x n) t_trend: datetime vector of computed trends
%   (n x m) trends: computed trends. m denotes the number of different time
%           series, i.e. E, N and U.
%   (1 x m) title_strings: title of every subplot for the m time series.
%   (u x 1) t_jump: datetime vector of jumps
%   (u x v) jump_categories_logical: array where row u = jump(u) and column
%           v = jump category(v). 1:=this jump (row) belongs to this 
%           category (column).
%   {1 x v} jump_categories_names: name of all categories, used in plot
%           legend.
%   (i x 1) jumpITRF: datetime vector of itrf jumps
%
%   OUTPUT:
%       Plot (Subplots when more than 1 time series component)

n_components = size(obs, 2); % number of trend components
jump_color_map = [...
    0, 204, 101; ...
    153, 51, 255; ...
    204, 0, 102; ...
    204, 204, 0; ...
    0, 204, 204; ...
    255, 128, 0
    ];

%% Create subplots for every coordinate component
for i = 1:n_components
    subplot(n_components, 1, i)
    
    % Plot Measurements
    pPts = plot(t_obs(~obs_outlier_logical{i}), obs(~obs_outlier_logical{i}, i), ...
        '.', 'markersize', 4, ...
        'Color', [0, 0.4470, 0.7410]);
    hold on

    % Plot Trend
    pTrend = plot(t_trend, trends(:, i), 'r');
    ylabel(sprintf('%s', obs_names{i}))
    xlim([min(t_obs) max(t_obs)])
    title(title_strings{i}, 'FontWeight','Normal');
    
    % Plot Outliers if there are any
    if nnz(obs_outlier_logical{i}) > 0
        % Outlier
        pOutl = plot(t_obs(obs_outlier_logical{i}), obs(obs_outlier_logical{i}, i), ...
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
    
    %% Set up plot legend
%     my_legend_elements = [pPts, pTrend, pITRF];
    my_legend_elements = [pPts, pTrend]; % w/o itrf
    % Update legend Items Outliers
    plotLogical = false(4 + length(jump_categories_names), 1);
    plotLogical(1:2) = true; % observations, trend,
%     plotLogical(1:3) = true; % observations, trend, itrf always in plot
    
    if nnz(obs_outlier_logical{i}) > 0 % if outliers present, add to legend
        plotLogical(4) = true;
        my_legend_elements = [my_legend_elements, pOutl];
    end
    
    plotElement = cell(length(jump_categories_names), 1); % plot elements storage
    cnt = 1; % counter for plot elements
    % Plot Jumps
    for n = 1:size(jump_categories_logical, 1) % All Jumps (=Rows)
        for c = 1:length(jump_categories_names) % All Categories (=Columns)
            if jump_categories_logical(n, c) % If this jumps fits to category
                
                plotElement{c} = plot([t_jump(n); t_jump(n)], [y1; y2], ... % plot jump
                    'color', jump_color_map(c, :)./255); % use color from map
                
                if ~plotLogical(4 + c) % if category of this jump is still FALSE ...
                    % this has do be done only ONCE per category, else
                    % there will be duplicate legend entries
                    my_legend_elements = [my_legend_elements, plotElement{c}]; % add jump category to legend
                end
                plotLogical(4 + c) = true; % set category of this jump to TRUE
            end
        end
    end
    all_plot_legend_names = [{'Observation', 'Trend', 'new ITRF', 'Outlier'}, jump_categories_names];
    my_plot_legend_names = all_plot_legend_names(plotLogical);
    
    % Add Legend
    legend(my_legend_elements, my_plot_legend_names, 'Location', 'eastoutside')
    hold off
end
end