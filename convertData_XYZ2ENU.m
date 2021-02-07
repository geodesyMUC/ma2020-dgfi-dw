clear variables
close all
addpath('myfunctions') % Add Function Storage to PATH

output_folder = 'data_psdENU';
% output_folder = 'data_psdXYZ'; % when processing XYZ (uncomment "XYZ data conversion") 

input_transl_file = 'src/PSD_GNSS';
input_folder = 'raw_data/Post-seismic-deformation/Daten/GNSS/';

myFolder_content.x = dir([input_folder, '*.x']); % station x
myFolder_content.y = dir([input_folder, '*.y']); % station y
myFolder_content.z = dir([input_folder, '*.z']); % station z

% other parameters
spheroid = referenceEllipsoid('GRS 80'); % itrf ellipsoid
t_origin = datetime('2000-01-01', 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');

% load all stations
for i = 1:length(myFolder_content.x) % y or z can also be used
    
    current_station_name = strrep(myFolder_content.x(i).name, '.x', ''); % name of station
%     current_station_nameSHRT = current_station_name(1:end-3); % cut last 3 chars
%     current_station_code = match_DOME_to_stationcode(current_station_nameSHRT, input_transl_file);
    % file paths
    current_station_path.x =  fullfile(myFolder_content.x(i).folder, myFolder_content.x(i).name);
    current_station_path.y =  fullfile(myFolder_content.y(i).folder, myFolder_content.y(i).name);
    current_station_path.z =  fullfile(myFolder_content.z(i).folder, myFolder_content.z(i).name);
    % observation data
    current_station_data.x = dlmread(current_station_path.x);
    current_station_data.y = dlmread(current_station_path.y);
    current_station_data.z = dlmread(current_station_path.z);
    
    % set up table for TSApart2A
    t_days = days(current_station_data.x(:, 1));
    % convert epochs to datetime
    d_t = t_origin + t_days; % [s] second intervals
    % calculate second intervals
    t = seconds(t_days - repmat(t_days(1), length(t_days), 1));
    x = current_station_data.x(:, 2); % [m]
    y = current_station_data.y(:, 2); % [m]
    z = current_station_data.z(:, 2); % [m]
    
    % convert xyz to lat lon h
    [lat,lon,h] = ecef2geodetic(spheroid, x, y, z, 'degrees');
    
    % convert lat lon h to enu
    lat0 = lat(1);  % origin latitude = first observation
    lon0 = lon(1);  % origin longitude = first observation
    h0 = h(1);      % origin height = first observation
    
    [E, N, U] = geodetic2enu(...
        lat,... % Latitude
        lon,... % Longitude
        h, ... % Ell. Height
        lat0, ... % Origin Latitude (use 1. obs)
        lon0, ... % Origin Longitude (use 1. obs)
        h0, ... % Origin Height (use 1. obs)
        spheroid, 'degrees'); % other parameters
    
    % load as table, convert enu to [mm]
    output_table = table(d_t, t, E .* 1e3, N .* 1e3, U .* 1e3, ...
    'VariableNames', ...
        {'date', 't', 'E', 'N', 'U'});

    % XYZ data conversion
    % the following code converts the XYZ time series in a format which can
    % be read by the main function without any conversion
    % the coordinates will be reduced by the first time series entry and converted to millimeters
%     output_table = table(d_t, t, (x-x(1)) .* 1e3, (y-y(1)) .* 1e3, (z-z(1)) .* 1e3, ...
%     'VariableNames', ...
%         {'date', 't', 'X', 'Y', 'Z'});
    
    currStation.Station = current_station_name;
    currStation.Data = output_table;
    
    % save output table
    save(fullfile(output_folder, [current_station_name, '.mat']), 'currStation');
    
    % create plot, save it
    figTS = figure('visible','off');
    VisualizeTS_ENU(currStation.Data, currStation.Station);
    saveas(figTS, fullfile(output_folder, ...
        [current_station_name, '-TSvis.png'])); % Save figure as image file
    
    % Set CreateFcn callback
    set(figTS, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    saveas(figTS, fullfile(output_folder, ...
        [current_station_name, '-TSvis.fig'])); % Save figure as image file
    close(figTS)
end

