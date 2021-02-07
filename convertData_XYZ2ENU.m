%CONVERTDATA_XYZ2ENU 
%    Converts formatted txt files with XYZ coordinates to .mat file which
%    are used as input for the main script. XYZ and ENU files and 
%    preliminary plots are created in the process.
%    This is the data pre-processing step.
%
%    Output is written to the specified directory. Note that ENU & XYZ will
%    have the same filename, but will be saved in different folders in
%    the output target directory.
%
%    (!) This script requires the MATLAP Mapping Toolbox.
%
% -------------------------------------------------------------------------
% Author: David Wallinger, DGFI. david.wallinger@live.at
% Master Thesis "Approximation of Non-Linear Post-Seismic Station 
% Motions in the Context of Geodetic Reference Frames" (2020)

clear variables
close all

addpath('myfunctions') % Add Function Storage to MATLAB path

% -------------------------------------------------------------------------
% Input paths
% Directory with x, y, and z coordinate files
input_folder = 'raw_data/Post-seismic-deformation/Daten/GNSS/';

% -------------------------------------------------------------------------
% Output directory
output_folder = 'data_psd';

% -------------------------------------------------------------------------
% Other parameters
spheroid = referenceEllipsoid('GRS 80'); % itrf ellipsoid
% Origin datetime
t_origin = datetime('2000-01-01', 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');

%% Process files
% Get all filepaths
myFolder_content.x = dir([input_folder, '*.x']); % station x
myFolder_content.y = dir([input_folder, '*.y']); % station y
myFolder_content.z = dir([input_folder, '*.z']); % station z

% load ALL stations
for i = 1:length(myFolder_content.x) % y or z can also be used
    % Get name of station from filename
    current_station_name = strrep(myFolder_content.x(i).name, '.x', '');
    currStation.Station = current_station_name;
    
    % ALternative (Deprecated)
    % input_transl_file = 'src/PSD_GNSS';
    % current_station_nameSHRT = current_station_name(1:end-3); % cut last 3 chars
    
    % Complete file path of current file
    current_station_path.x =  fullfile(myFolder_content.x(i).folder, myFolder_content.x(i).name);
    current_station_path.y =  fullfile(myFolder_content.y(i).folder, myFolder_content.y(i).name);
    current_station_path.z =  fullfile(myFolder_content.z(i).folder, myFolder_content.z(i).name);
    % Read the delimited data
    current_station_data.x = dlmread(current_station_path.x);
    current_station_data.y = dlmread(current_station_path.y);
    current_station_data.z = dlmread(current_station_path.z);
    
    %% Set up table
    % Convert time in days since origin to datetime
    t_days = days(current_station_data.x(:, 1));
    d_t = t_origin + t_days; % [s] second intervals
    
    % Convert time in seconds
    t = seconds(t_days - repmat(t_days(1), length(t_days), 1));
    
    % Coordinate conversion
    x = current_station_data.x(:, 2); % [m]
    y = current_station_data.y(:, 2); % [m]
    z = current_station_data.z(:, 2); % [m]
    
    %% Get ENU coordinates
    % Convert XYZ to geodetic coordinates lat, lon, h (LLH)
    [lat,lon,h] = ecef2geodetic(spheroid, x, y, z, 'degrees');
    
    % Convert LLH to ENU relative to FIRST observation
    lat0 = lat(1);  % origin latitude = first observation
    lon0 = lon(1);  % origin longitude = first observation
    h0 = h(1);      % origin height = first observation
    scaleF = 1e3;   % [mm]    
    
    [E, N, U] = geodetic2enu(...
        lat,...               % Latitude
        lon,...               % Longitude
        h, ...                % Ell. Height
        lat0, ...             % Origin Latitude (use 1. obs)
        lon0, ...             % Origin Longitude (use 1. obs)
        h0, ...               % Origin Height (use 1. obs)
        spheroid, 'degrees'); % other parameters
    
    % Load in table, convert ENU to [mm]
    output_tableENU = table(d_t, t, E .* scaleF, N .* scaleF, U .* scaleF, ...
    'VariableNames', ...
        {'date', 't', 'E', 'N', 'U'});
    currStation.Data = output_tableENU;
    
    
    
    % create plot
    figTS = figure('visible','off');
    VisualizeTS_ENU(currStation.Data, currStation.Station);
    % save output
    if ~exist(fullfile(output_folder, 'ENU'), 'dir')
        mkdir(fullfile(output_folder, 'ENU'))
    end
    if ~exist(fullfile(output_folder, 'ENU_plots'), 'dir')
        mkdir(fullfile(output_folder, 'ENU_plots'))
    end
    save(fullfile(output_folder, 'ENU', [current_station_name, '', '.mat']), 'currStation');
    saveas(figTS, fullfile(output_folder, 'ENU_plots',...
        [current_station_name, '-ENU-TSvis.png'])); % Save figure as image file
    close(figTS)
    
    %% Get XYZ coordinates
    % Relative to first observation in MM
    output_tableXYZ = table(d_t, t, ...
        (x-x(1)) .* scaleF, (y-y(1)) .* scaleF, (z-z(1)) .* scaleF, ...
        'VariableNames', ...
        {'date', 't', 'X', 'Y', 'Z'});
    currStation.Data = output_tableXYZ;
    % Create plot
    figTS = figure('visible','off');
    VisualizeTS_ENU(currStation.Data, currStation.Station);
    % Save output
    if ~exist(fullfile(output_folder, 'XYZ'), 'dir')
        mkdir(fullfile(output_folder, 'XYZ'))
    end
    if ~exist(fullfile(output_folder, 'XYZ_plots'), 'dir')
        mkdir(fullfile(output_folder, 'XYZ_plots'))
    end
    save(fullfile(output_folder, 'XYZ', [current_station_name, '', '.mat']), 'currStation');
    saveas(figTS, fullfile(output_folder, 'XYZ_plots',  ...
        [current_station_name, '-XYZ-TSvis.png'])); % Save figure as image fil
    close(figTS)
end

