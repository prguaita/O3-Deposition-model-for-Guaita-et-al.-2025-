%% Main Script for Ozone Deposition and Plant-Atmosphere Interaction Model
%
% This script simulates ozone deposition and stomatal conductance dynamics
% over a specified period, using meteorological, soil, and plant physiological
% inputs. The model accounts for light and accumulation periods, integrates
% Jarvis functions for stomatal responses, and calculates ozone effects on
% wheat canopy at multiple heights.
%
% Key features:
%   - Accumulation of stomatal conductance and ozone concentration metrics
%     during light hours for POD accumulation periods (PAC).
%   - Dynamic interpolation of input maps over time.
%   - Calculation of stomatal conductance fractions for sunlit and shaded leaves
%     with different parameterizations (winter, Mediterranean).
%   - Soil water content and plant-available water modeled or read from data.
%   - Annual saving of plant ozone dose (POD) maps and related diagnostic variables.
%   - Support for multiple CMIP6 experiments configurations
%
% Inputs (loaded or defined externally):
%   - Meteorological time series (temperature, VPD, radiation)
%   - Soil moisture indices and hydraulic parameters
%   - Plant physiological parameters (LAI, thermal time, conductance params)
%   - Ozone concentration data
%
% Outputs:
%   - Maps (or Time series) of stomatal conductance metrics
%   - Annual POD maps
%   - Diagnostic tables summarizing key variables per timestep
%   - Saved data files for further analysis or use in subsequent model runs
%
% Usage:
%   Configure experiment parameters in 'config' struct.
%   Run the script to simulate for the defined period and save results.
%
% Author: PR Guaita
% Date (first version for global study): July 2023
%
% Notes:
%   - Ensure supporting functions and data files are available in the MATLAB path.
%   - Adjust 'folder_save', 'folder_main', and other path variables as needed.
%   - Use config.table_out_option to enable or disable diagnostic output tables (only for single-node runs!).
%
% Documentation for the non-global version:
%   - Guaita, P. R., Marzuoli, R., and Gerosa, G. A.: A regional scale flux-based O3 risk assessment for winter wheat 
%   in northern Italy, and effects of different spatio-temporal resolutions, Environmental Pollution, 333, 121860, 
%   https://doi.org/10.1016/j.envpol.2023.121860, 2023.
%

clear;  % Clear workspace

%% Load configuration and define paths

config       = ini2struct('config_POD_wheat.ini');   % Read config file
file_suffix  = config.file_suffix;
disp(['-----------------------' file_suffix])

folder_main      = config.folder_main;
model_name       = config.model_name;
addpath(genpath(folder_main))  % Add main folder and subfolders to MATLAB path

folder_result    = fullfile(config.folder_result, model_name, file_suffix);
folder_data      = fullfile(config.folder_data, model_name, file_suffix);
folder_save      = fullfile(folder_main, config.folder_result, model_name, ...
                    config.plant_species, file_suffix, config.experiment, ...
                    ['hem_' config.filter_hem]);

% Create output folder if it doesn't exist
if ~isfolder(folder_save)
    mkdir(folder_save);
end

%% Load stomatal conductance parameters

par_g_stom_winter = eval(['config.g_par_' config.dep.conductance_type]);
par_g_stom_spring = eval(['config.g_par_' config.dep.conductance_type '_spring']);
par_g_stom_medite = eval(['config.g_par_' config.dep.conductance_type '_mediterranean']);

%% Simulation period and days per year

yr_start    = config.yr_start;
yr_end      = config.yr_end;
year_array  = yr_start:yr_end;
n_day_year  = ones(size(year_array));   % Initialize array

switch config.input_model.opt_calendar
    case 'fix'
        n_day_year = config.input_model.n_day_year * n_day_year;
    case 'leap'
        n_day_year(leapyear(year_array))     = 366;
        n_day_year(~leapyear(year_array))    = 365;
end

%% Accumulation periods in degree days from sowing

A_start            = config.plant_TT.TT_A_start;
A_end              = config.plant_TT.TT_A_end;
A_start_spring     = config.plant_TT_spring.TT_A_start;
A_end_spring       = config.plant_TT_spring.TT_A_end;
A_start_medit      = config.plant_TT_mediterranean.TT_A_start;
A_end_medit        = config.plant_TT_mediterranean.TT_A_end;

%% Maximum conductances (for water and ozone)

g_max_W_winter     = par_g_stom_winter.g_max_W;
g_max_W_medite     = par_g_stom_medite.g_max_W;

% Scale by diffusivity ratio for ozone
g_max_O3_winter    = g_max_W_winter * config.const.D01_O3 / config.const.D01_W;
g_max_O3_medite    = g_max_W_medite * config.const.D01_O3 / config.const.D01_W;

%% Load static maps

tmp_suffix = [];
if ~isempty(file_suffix)
    tmp_suffix = ['_' file_suffix];
end

% Latitude and longitude grids
load(fullfile(folder_main, folder_data, 'aux_data', ...
     [model_name '_static_horgrid' tmp_suffix '.mat']));

natlat = lat;
natlon = lon;

%% Filter latitude based on hemisphere

switch config.filter_hem
    case 'N'
        flag_lat = lat > 0;
        lat_hem  = lat(flag_lat);
    case 'S'
        flag_lat = lat <= 0;
        lat_hem  = lat(flag_lat);
    case 'p'
        flag_lat = true(size(lat));
        lat_hem  = lat;
        timezone_hr = point_table.local_hr(i_point); % Note: placeholder
end

%% Load sowing dates

load(fullfile(folder_main, folder_data, 'aux_data', ...
    [model_name '_' config.experiment '_sowing_dates_wheat' tmp_suffix '.mat']));

switch config.experiment
    case 'historical'
        sow_map = sow_map(:, :, (yr_start - 1999):(yr_end - 1999));
    case 'ssp370pdSST'
        sow_map = sow_map(:, :, (yr_start - 2089):(yr_end - 2089));
    otherwise
        sow_map = sow_map(:, :, (yr_start - 2014):(yr_end - 2014));
end

%% Load climate classification and spring wheat flag

load(fullfile(folder_main, folder_data, 'aux_data', ...
    [model_name '_static_wheat_climate' tmp_suffix '.mat']));

n_row = length(lon);
n_col = length(lat_hem);

%% Build model time arrays

n_hours_year = n_day_year * 24;
timestep     = config.timestep;  % timestep in hours

% Build array of fractional days across all years
time_tgt = [];
for i_year = 1:length(year_array)
    time_tmp = 1 + (0:timestep:n_hours_year(i_year)) / 24;
    time_tmp(end) = [];  % Remove last to keep consistent length
    time_tgt = [time_tgt, time_tmp];
end

%% Build cumulative days and datetime arrays

day_year_array_greg_cumul = [0, cumsum(365 + leapyear(year_array))];
day_year_array_greg_cumul(end) = [];

day_year_array_greg = [];
for i_year = 1:length(n_hours_year)
    n_steps = n_hours_year(i_year) / timestep;
    day_year_array_greg = [day_year_array_greg, ...
        repmat(day_year_array_greg_cumul(i_year), 1, n_steps)];
end

% Full datetime array
date_array = datetime(yr_start,1,1,0,0,0) + day_year_array_greg + time_tgt - 1;

%% Build CMIP-style time coordinates

[time_start, time_end] = timegreg_2_timeCMIP( ...
    datetime(yr_start,1,1), datetime(yr_end+1,1,1), datetime(yr_start,1,1), model_name);

[time_start_ref, ~] = timegreg_2_timeCMIP( ...
    datetime(yr_start,1,1), datetime(yr_end+1,1,1), datetime(1850,1,1), model_name);

% Array of timesteps in days
time_tgt_days = time_start : (timestep/24) : time_end;
time_tgt_days(end) = [];

%% Apply hemisphere filter to data maps

sow_map(:, ~flag_lat, :)          = [];
whclim_map(:, ~flag_lat, :)       = [];
spring_wheat_map(:, ~flag_lat, :) = [];

%% Load meteorological data and PLANT PHENOLOGY module
% This section loads various input datasets required for the model,
% including meteorological variables and plant phenology inputs.
% Because of the large data volume and processing, this step can take
% around two hours to complete depending on hardware.

tic  % Start timing

%% Load ozone geometric height map (z_m_O3)

% Identify the file containing ozone height map data for the experiment
file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_z_*']));
load(fullfile(file_path.folder, file_path.name));

% Extract the geometric height and adjust time reference to simulation start
z_m_O3 = geom_z;
z_m_O3_time = time_serie - time_start_ref + 1;

if numel(z_m_O3_time) > 1
    % If multiple time points exist, filter times within model simulation period
    flag_time_z = time_start <= z_m_O3_time & z_m_O3_time <= time_end;
    
    % Remove times and corresponding height maps outside this period
    z_m_O3_time(~flag_time_z) = [];
    z_m_O3(:, :, ~flag_time_z) = [];

    % Prepare grids for interpolation of height data onto native model grid filtered by hemisphere
    [geomgrid_lon, geomgrid_lat, geomgrid_time] = ndgrid(lon, lat, z_m_O3_time);
    geomint = griddedInterpolant(geomgrid_lon, geomgrid_lat, geomgrid_time, z_m_O3, 'nearest', 'nearest');
    [natgrid_lon, natgrid_lat, natgeomgrid_time] = ndgrid(lon, lat_hem, z_m_O3_time);
    
    % Interpolate height data to native grid points
    z_m_O3 = geomint(natgrid_lon, natgrid_lat, natgeomgrid_time);
else
    % If height data is constant in time (only one map), just regrid spatially
    [geomgrid_lon, geomgrid_lat] = ndgrid(lon, lat);
    geomint = griddedInterpolant(geomgrid_lon, geomgrid_lat, z_m_O3, 'nearest', 'nearest');
    [natgrid_lon, natgrid_lat] = ndgrid(lon, lat_hem);
    z_m_O3 = geomint(natgrid_lon, natgrid_lat);
end
disp('Loaded geometric height (geom_z)')

%% Load temperature data (in Kelvin) and convert to Celsius

% Identify temperature file for experiment and period
file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_tas_' int2str(yr_start) '-' int2str(yr_end) '*']));
load(fullfile(file_path.folder, file_path.name));

% Filter spatial data for the hemisphere of interest
mat_data(:, ~flag_lat, :) = [];

% Convert temperature from Kelvin to Celsius for biological relevance
tas_C = mat_data - 273.2;

% Adjust time reference relative to dataset start date
tas_time = time_serie - floor(min(time_serie)) + 1;

disp(['Loaded temperature data: ' file_path.name])
clear mat_data

%% Calculate thermal time maps needed for plant phenology

% Create a flag matrix marking when to start accumulating thermal time based on sowing date
flag_TT_init = (cumsum(ones(n_row, n_col, n_day_year(1)), 3) / n_day_year(1)) > sow_map(:, :, 1);

% Find number of timesteps in a day, assuming uniform time spacing
n_timestep_day = find(floor(tas_time) > 1, 1) - 1;

% Calculate daily average temperature, clipped at plant minimum thermal threshold,
% and consider only days after sowing
tas_day = max(config.plant_TT.TT_min, ...
    mean(reshape(tas_C(:, :, 1:(n_hours_year/24 * n_timestep_day)), ...
    size(tas_C, 1), size(tas_C, 2), [], n_timestep_day), 4, 'omitnan'));

% Apply the sowing date mask to ignore thermal time before sowing
tas_day = flag_TT_init .* tas_day;

% Initialize thermal time accumulation map by summing daily values
TT_map_init = sum(tas_day, 3);

% Calculate accumulated thermal time over the full period using custom function
TT_data = thermal_time_v2(tas_C, tas_time, config, sow_map, TT_map_init);

%% Run PLANT PHENOLOGY module

% Using accumulated thermal time and climate data, compute plant state variables:
% LAI = Leaf Area Index, SAI = Surface Area Index, RD = Root Depth,
% Alb = Albedo, h_c = Canopy height, d = Diameter, z0M and z0H = roughness lengths
[LAI_map, SAI_map, RD_map, Alb_map, h_c_map, ...
    d_map, z0M_map, z0H_map] = module_plant_phenology_v1(TT_data, whclim_map, config);

%% Load ozone concentration [ppb]

file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_sfo3_' int2str(yr_start) '-' int2str(yr_end) '*']));
load(fullfile(file_path.folder, file_path.name));

% Filter for hemisphere and store ozone data and adjusted time vector
mat_data(:, ~flag_lat, :) = [];
sfo3_ppb = mat_data;
sfo3_time = time_serie - floor(min(time_serie)) + 1;

disp(['Loaded ozone concentration: ' file_path.name])
clear mat_data

%% Load surface pressure [Pa] and convert to kilopascals (kPa)

file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_ps_' int2str(yr_start) '-' int2str(yr_end) '*']));
load(fullfile(file_path.folder, file_path.name));
mat_data(:, ~flag_lat, :) = [];

ps_kPa = round(mat_data / 1000, 2);  % Convert Pa to kPa with rounding
ps_time = time_serie - floor(min(time_serie)) + 1;

disp(['Loaded surface pressure: ' file_path.name])
clear mat_data

%% Load specific humidity (dimensionless) and convert to vapor pressure [kPa]

file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_huss_' int2str(yr_start) '-' int2str(yr_end) '*']));
load(fullfile(file_path.folder, file_path.name));
mat_data(:, ~flag_lat, :) = [];

% Find index of closest matching surface pressure time for each specific humidity timestep
closestElements = zeros(size(time_serie));
for i_timestep = 1:length(time_serie)
    [~, idx] = min(abs(ps_time - time_serie(i_timestep)));
    closestElements(i_timestep) = idx;
end

% Convert specific humidity q to vapor pressure e using formula: e = p * q / 0.622
switch ndims(mat_data)
    case 2
        e_kPa = ps_kPa(closestElements) .* mat_data / 0.622;
    case 3
        e_kPa = ps_kPa(:, :, closestElements) .* mat_data / 0.622;
end

e_kPa_time = time_serie - floor(min(time_serie)) + 1;

disp(['Loaded specific humidity: ' file_path.name])
clear mat_data

%% Load precipitation data [mm/day] and convert to mm per timestep

file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_pr_' int2str(yr_start) '-' int2str(yr_end) '*']));
load(fullfile(file_path.folder, file_path.name));
mat_data(:, ~flag_lat, :) = [];

% Convert precipitation to rate per model timestep (hours to days conversion)
pr_mm = mat_data / (24 / timestep);
pr_time = time_serie - floor(min(time_serie)) + 1;

disp(['Loaded precipitation: ' file_path.name])
clear mat_data

%% Load surface wind speed [m/s]

file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_sfcWind_' int2str(yr_start) '-' int2str(yr_end) '*']));
load(fullfile(file_path.folder, file_path.name));
mat_data(:, ~flag_lat, :) = [];

sfcWind_ms = mat_data;
sfcWind_time = time_serie - floor(min(time_serie)) + 1;

disp(['Loaded wind speed: ' file_path.name])
clear mat_data

%% Load global shortwave radiation [W/m²]

file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_rsds_' int2str(yr_start) '-' int2str(yr_end) '*']));
load(fullfile(file_path.folder, file_path.name));
mat_data(:, ~flag_lat, :) = [];

rsds_Wm2 = mat_data;
rsds_time = time_serie - floor(min(time_serie)) + 1;

disp(['Loaded shortwave radiation: ' file_path.name])
clear mat_data

%% Load atmospheric CO2 concentration [ppm]

if ~strcmp(config.experiment, 'ssp370pdSST')
    % For standard experiments, load experiment-specific CO2 data
    file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_co2_' int2str(yr_start) '-' int2str(yr_end) '*']));
    load(fullfile(file_path.folder, file_path.name));
    co2_ppm = mat_data;
    co2_time = time_serie - floor(min(time_serie)) + 1;
else
    % For SSP370pdSST, load historical CO2 data as fallback
    file_path = dir(fullfile(folder_main, folder_data, ['*historical_co2_*']));
    load(fullfile(file_path.folder, file_path.name));
    co2_ppm = mat_data(1:132);  % Limit to subset if necessary
    co2_time = time_serie - floor(min(time_serie)) + 1;
end

disp(['Loaded CO2 concentration: ' file_path.name])
clear mat_data

%% Clear large temporary variables to free memory

clear geom_z mean_z prct_5 prct_50 prct_95 lat lon time_serie natgeomgrid_time ...
    geomgrid_lon geomgrid_lat geomgrid_time geomint mat_data

%% Load optional surface energy balance variables: net radiation and sensible heat flux

% Net radiation (Rn) loading and filtering
if config.input_model.Rn
    file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_rnds_' int2str(yr_start) '-' int2str(yr_end) '*']));
    load(fullfile(file_path.folder, file_path.name));
    mat_data(:, ~flag_lat, :) = [];
    Rn_Wm2 = mat_data;
    Rn_time = time_serie - floor(min(time_serie)) + 1;
else
    % If not provided, initialize as NaNs for the entire time series
    Rn_Wm2 = squeeze(nan([size(sow_map, 1), size(sow_map, 2), length(time_tgt_days)]));
    Rn_time = time_tgt_days;
end
disp(['Loaded net radiation: ' file_path.name])
clear mat_data

% Sensible heat flux (H) loading and filtering
if config.input_model.H
    file_path = dir(fullfile(folder_main, folder_data, ['*' config.experiment '_hfss_' int2str(yr_start) '-' int2str(yr_end) '*']));
    load(fullfile(file_path.folder, file_path.name));
    mat_data(:, ~flag_lat, :) = [];
    H_Wm2 = mat_data;
    H_time = time_serie - floor(min(time_serie)) + 1;
else
    % If not provided, initialize as NaNs
    H_Wm2 = squeeze(nan([size(sow_map, 1), size(sow_map, 2), length(time_tgt_days)]));
    H_time = time_tgt_days;
end
disp(['Loaded sensible heat flux: ' file_path.name])
clear mat_data

%% Report total loading time

disp(['Loading input variables completed in ' num2str(toc/60, '%.2f') ' minutes.'])

%% Solar elevation [degree] calculation ------------------------------
% This section calculates the solar elevation angle for each grid point
% and each timestep over the simulation period. It is a computationally 
% intensive step and typically takes 20-30 minutes to complete, depending 
% on hardware.

% Initialize the solar elevation matrix B depending on calendar type
switch config.input_model.opt_calendar
    case 'fix'
        % For fixed calendar: define matrix dimensions based on number of longitudes,
        % latitudes, and total number of timesteps (hours/timestep * days in year)
        B = zeros(length(natlon), length(lat_hem), (24 / timestep) * config.input_model.n_day_year);
    case 'leap'
        % For leap year calendar: use 365 days for total timesteps
        B = zeros(length(natlon), length(lat_hem), (24 / timestep) * 365);
end

% Adjust longitude values to match range expected by SolarAzEl function
% SolarAzEl expects longitude in range [-180, 180], so subtract 360 for values >= 180
lon_B = natlon;
lon_B(lon_B >= 180) = lon_B(lon_B >= 180) - 360;

% Select a subset of latitudes and longitudes to reduce computation during initial calculation
flag_lat_B = (1:7:length(lat_hem))';  % every 7th latitude index
flag_lon_B = (1:10:length(lon_B))';   % every 10th longitude index

% Ensure the last latitude and longitude index are included in flags for full coverage
if lat_hem(flag_lat_B(end)) ~= lat_hem(end)
    flag_lat_B(end + 1) = length(lat_hem);
end
if lon_B(flag_lon_B(end)) ~= lon_B(end)
    flag_lon_B(end + 1) = length(lon_B);
end

% Create mesh grids for all native longitude and latitude points
[natgrid_lon, natgrid_lat] = ndgrid(natlon, lat_hem);

% Create mesh grids for the reduced subset of points used for initial solar elevation calculation
[gridB_lon, gridB_lat] = ndgrid(natlon(flag_lon_B), lat_hem(flag_lat_B));

% Initialize a temporary matrix to store solar elevation values at the reduced grid points
tmp_B = zeros(length(flag_lon_B), length(flag_lat_B));

tic  % Start timing the solar elevation calculations

% Loop over each timestep in the simulation period
for i_timestep = 1:size(B, 3)
    % Loop over the reduced longitude subset
    for i_lon = 1:length(flag_lon_B)
        % Loop over the reduced latitude subset
        for i_lat = 1:length(flag_lat_B)
            % Extract the longitude and latitude for this grid point
            tmp_lon = lon_B(flag_lon_B(i_lon));
            tmp_lat = lat_hem(flag_lat_B(i_lat));
            
            % Calculate solar azimuth and elevation angles for this date and location
            % SolarAzEl returns azimuth and elevation; we only keep elevation (second output)
            switch ndims(tas_C)
                case 2
                    [~, tmp_B] = SolarAzEl(date_array(i_timestep), tmp_lat, tmp_lon, 0);
                case 3
                    [~, tmp_B(i_lon, i_lat)] = SolarAzEl(date_array(i_timestep), tmp_lat, tmp_lon, 0);
            end
        end
    end
    
    % After calculating for reduced grid, interpolate to full native grid
    switch ndims(tas_C)
        case 2
            % If tas_C is 2D, just assign the single value (scalar)
            B(i_timestep) = tmp_B(i_timestep);
        case 3
            % If tas_C is 3D, interpolate from reduced grid (tmp_B) to full grid (natgrid_lon, natgrid_lat)
            B_interp = griddedInterpolant(gridB_lon, gridB_lat, tmp_B, 'spline', 'spline');
            B(:, :, i_timestep) = B_interp(natgrid_lon, natgrid_lat);
    end
end

% Limit solar elevation values to zero as minimum, since negative solar elevation
% means sun is below horizon and does not contribute to solar radiation
B(B < 0) = 0;

disp(['Solar elevation calculation completed in ' num2str(toc/60) ' minutes.'])

%% Initializing SOIL module ---------------------------------------------

% Load soil-related parameters necessary for modeling soil water dynamics.
% These parameters include Field Capacity (FC), Plant Available Water (Lm3),
% and Wilting Point (WP). They define key soil moisture thresholds affecting plant water availability.

% Load Field Capacity (FC) from auxiliary data file
% FC represents the amount of water soil can hold against gravity (in volumetric units).
load(fullfile(folder_main, folder_data, 'aux_data', 'Field_capacity_mean.mat'));
FC = var_data;  % Assign loaded variable to FC

% Load Plant Available Water (Lm3) from auxiliary data file
% Lm3 here is converted from the loaded variable (assumed in meters) to liters per cubic meter
load(fullfile(folder_main, folder_data, 'aux_data', 'Plant_available_water_mean.mat'));
Lm3 = var_data * 1000;  % Convert to liters per cubic meter (L/m³)

% Calculate Wilting Point (WP) as the difference between Field Capacity and
% Plant Available Water, representing soil moisture level below which plants cannot extract water.
WP = FC - Lm3 / 1000;  % Convert Lm3 back to volumetric units before subtraction

% Filter data based on longitude selection mask 'flag_lat'
% This removes data columns for longitudes that are not in the selected latitude range.
WP(:, not(flag_lat)) = [];
FC(:, not(flag_lat)) = [];
Lm3(:, not(flag_lat)) = [];

%% Initializing variables ---------------------------------------------------

% Pre-allocate arrays for model state variables and diagnostic outputs.
% Most arrays are sized according to the spatial dimensions of 'sow_map' (e.g., lat x lon grid).
% This improves performance by avoiding dynamic resizing during simulation.

% Initialize zero matrices for various physiological, environmental, and soil variables
sumVPD             = zeros(size(sow_map,[1,2]));  % Accumulated Vapor Pressure Deficit
g_s_frac_sun_prev  = zeros(size(sow_map,[1,2]));  % Previous stomatal conductance fraction (sunlit leaves)
g_s_frac_shad_prev = zeros(size(sow_map,[1,2]));  % Previous stomatal conductance fraction (shaded leaves)
light_flag_prev    = zeros(size(sow_map,[1,2]));  % Previous light condition flag

AWC_prev           = Lm3;                           % Available Water Content initialized from soil data
RD                 = zeros(size(sow_map,[1,2]));  % Root depth
RD_prev            = zeros(size(sow_map,[1,2]));  % Previous root depth
LAI_prev           = zeros(size(sow_map,[1,2]));  % Previous Leaf Area Index

W_in_prev          = zeros(size(sow_map,[1,2]));  % Previous soil water input
S_c                = zeros(size(sow_map,[1,2]));  % Soil moisture stress coefficient

% Evaporation and transpiration storages (previous timestep)
evap_W_prev        = zeros(size(sow_map,[1,2]));
tras_W_prev        = zeros(size(sow_map,[1,2]));

% Physiological Ozone Damage (POD) accumulators and previous values
POD0               = zeros(size(sow_map,[1,2]));  % POD baseline
PODY               = zeros(size(sow_map,[1,2]));  % POD yield
POD0_prev          = zeros(size(sow_map,[1,2]));  % Previous POD0
PODY_prev          = zeros(size(sow_map,[1,2]));  % Previous PODY

% Final accumulations of POD metrics
PODY_final         = zeros(size(sow_map,[1,2]));
POD0_final         = zeros(size(sow_map,[1,2]));

% Counters and flags for Photochemical Activity Cycles (PAC)
n_timestep_PAClight = zeros(size(sow_map,[1,2]));

% Outputs from Jarvis stomatal conductance model components (temperature, soil moisture, VPD)
jarvis_out_temp    = zeros(size(sow_map,[1,2]));
jarvis_out_soil    = zeros(size(sow_map,[1,2]));
jarvis_out_VPD     = zeros(size(sow_map,[1,2]));

% Surface ozone concentration variables (ppb = parts per billion)
sfo3_ppb_out       = zeros(size(sow_map,[1,2]));
sfo3_diff_ppb      = zeros(size(sow_map,[1,2]));

% Start day for PAC calculation (initialized as NaN)
start_day_PAC      = nan * ones(size(sow_map,[1,2]));

% Additional environmental output variables (temperature, vapor pressure deficit, plant available water)
T_out              = zeros(size(sow_map,[1,2]));
VPD_out            = zeros(size(sow_map,[1,2]));
PAW_out            = zeros(size(sow_map,[1,2]));

% Final aggregated outputs from Jarvis model and ozone variables
jarvis_out_temp_final     = zeros(size(sow_map,[1,2]));
jarvis_out_soil_final     = zeros(size(sow_map,[1,2]));
jarvis_out_VPD_final      = zeros(size(sow_map,[1,2]));
sfo3_ppb_out_final       = zeros(size(sow_map,[1,2]));
sfo3_diff_ppb_final      = zeros(size(sow_map,[1,2]));
T_out_final              = zeros(size(sow_map,[1,2]));
VPD_out_final            = zeros(size(sow_map,[1,2]));
PAW_out_final            = zeros(size(sow_map,[1,2]));
n_day_PAC_final          = zeros(size(sow_map,[1,2]));
start_day_PAC_final      = nan * ones(size(sow_map,[1,2]));

% If running historical experiment, allocate large 3D arrays for ozone output
if strcmp(config.experiment, 'historical')
    o3_2m_output_ppb = single(nan(length(lon), length(lat_hem), length(time_tgt)));
    o3_3m_output_ppb = single(nan(length(lon), length(lat_hem), length(time_tgt)));
end

% For non-historical experiments, load previous simulation data for iterative calculations
if ~strcmp('historical', config.experiment)
    filename = fullfile(folder_main, folder_data, ...
        [model_name '_historical_prevData_' config.filter_hem '_2014.mat']);
    load(filename);
end

% Initialize diagnostic output table if requested in config
if config.table_out_option
    % Variable names to be stored in the diagnostic table for each timestep
    var_names_array = {'yr', 'doy', 'mth', 'hr', 'hr_local', ...
        'tas_tmp', 'sfcWind_tmp', 'rsds_tmp', 'ps_tmp', 'e_tmp', 'sfo3_tmp', 'sfo3_zm_ugm3', ...
        'pr_tmp', 'LE_tmp', 'Rn_tmp', 'H_tmp', 'B_tmp', 'VPD_kPa', ...
        'LAI', 'h_c', 'PAC_flag', 'TT', ...
        'MWHC', 'AWC', 'SWC_tmp', 'PAW_tmp', ...
        'light_flag', ...
        'g_s_frac_sun', 'jarvis_func.fphen', 'jarvis_func.ftemp', 'jarvis_func.fVPD', 'jarvis_func.fsoil', ...
        'u_star_hour', 'L_hour', 'u_hc', 'R_aH_hc_zmO3', 'R_bO3', 'R_surf_O3', 'R_aH_dz0m_zmT', 'R_inc', ...
        'r_s_act_O3', 'rc', 'rbO3', 'g_ext', 'O3_h_c_ppb', 'F_O3', 'PODY', 'POD0', ...
        'tras_W', 'E_soil', 'E_wet', 'S_c', 'W_in'};

    % Define all variables as double precision for table storage
    var_type_array = repmat({'double'}, 1, length(var_names_array));

    % Pre-allocate diagnostic table with appropriate size and variable types
    diag_table = table('Size', [length(time_tgt_days), length(var_names_array)], ...
                       'VariableNames', var_names_array, ...
                       'VariableTypes', var_type_array);
end

%% run the model

tic
for i_timestep = 1:length(time_tgt_days)

    %% Time step preprocessing 
    
    % Extract calendar variables from the current timestep's date
    yr  = year(date_array(i_timestep));                  % Year
    doy = day(date_array(i_timestep), 'dayofyear');      % Day of year
    mth = month(date_array(i_timestep));                 % Month
    hr  = hour(date_array(i_timestep));                   % Hour
    
    % Calculate day of the period (dop) relative to simulation start
    % Example: If calculations start in 2015, Jan 1 2015 corresponds to dop = 1
    dop = floor(time_tgt_days(i_timestep));
    
    % Interpolate meteorological and environmental variables for the current timestep
    % Nearest neighbor interpolation is used to minimize computational cost.
    tas_tmp     = interp_time(time_tgt_days(i_timestep), tas_time, tas_C, 'nearest');       % Air temperature (°C)
    sfcWind_tmp = interp_time(time_tgt_days(i_timestep), sfcWind_time, sfcWind_ms, 'nearest'); % Surface wind speed (m/s)
    rsds_tmp    = interp_time(time_tgt_days(i_timestep), rsds_time, rsds_Wm2, 'nearest');   % Downward shortwave radiation (W/m²)
    ps_tmp      = interp_time(time_tgt_days(i_timestep), ps_time, ps_kPa, 'nearest');       % Surface pressure (kPa)
    e_tmp       = interp_time(time_tgt_days(i_timestep), e_kPa_time, e_kPa, 'nearest');     % Vapor pressure (kPa)
    sfo3_tmp    = interp_time(time_tgt_days(i_timestep), sfo3_time, sfo3_ppb, 'nearest');   % Surface ozone (ppb)
    pr_tmp      = interp_time(time_tgt_days(i_timestep), pr_time, pr_mm, 'nearest');        % Precipitation (mm)
    Rn_tmp      = interp_time(time_tgt_days(i_timestep), Rn_time, Rn_Wm2, 'nearest');       % Net radiation (W/m²)
    H_tmp       = interp_time(time_tgt_days(i_timestep), H_time, H_Wm2, 'nearest');         % Sensible heat flux (W/m²)
    z_m_O3_tmp  = interp_time(time_tgt_days(i_timestep), z_m_O3_time, z_m_O3, 'nearest');   % Measurement height for O3 (m)
    co2_tmp     = interp_time(time_tgt_days(i_timestep), co2_time, co2_ppm, 'nearest');     % CO2 concentration (ppm)
    
    % Determine the base timestep index for climatological forcing variables
    switch config.input_model.opt_calendar
        case 'fix'
            % Fixed calendar: wrap around based on hours per year for the given year
            B_timestep = mod(i_timestep - 1, n_hours_year(yr - yr_start + 1)) + 1;
        case 'leap'
            % Leap year aware: calculate based on day of year capped at 365, and hourly timestep
            B_timestep = (min(365, day(date_array(i_timestep), 'dayofyear')) - 1) * (24 / timestep) ...
                         + mod(i_timestep - 1, 24 / timestep) + 1;
    end
    
    % Extract the relevant forcing map B at the current timestep
    switch ndims(tas_C)
        case 2
            B_tmp = B(B_timestep);          % 2D forcing
        case 3
            B_tmp = B(:, :, B_timestep);    % 3D forcing (e.g., lat x lon x time)
    end
    
    % Load daily vegetation and surface data depending on data dimensionality
    switch ndims(tas_C)
        case 2
            % 2D data maps indexed by day of period (dop)
            LAI  = LAI_map(dop);
            SAI  = SAI_map(dop);
            RD   = RD_map(dop);
            h_c  = h_c_map(dop);
            d    = d_map(dop);
            z_0m = z0M_map(dop);
            z_0H = z0H_map(dop);
            Alb  = Alb_map(dop);
    
            % Determine photochemical accumulation period flags based on temperature thresholds
            PAC_flag = (TT_data(dop) > A_start) & (TT_data(dop) < A_end);
            PAC_flag_spring = (TT_data(dop) > A_start_spring) & (TT_data(dop) < A_end_spring);
            PAC_flag_medit  = (TT_data(dop) > A_start_medit) & (TT_data(dop) < A_end_medit);
    
            % Adjust PAC flags for specific climate regions
            PAC_flag(whclim_map == 1.5)  = PAC_flag_spring(whclim_map == 1.5);
            PAC_flag(whclim_map == 0.5)  = PAC_flag_medit(whclim_map == 0.5);
            PAC_flag(whclim_map == -0.5) = PAC_flag_medit(whclim_map == -0.5);
    
            TT = TT_data(dop);
    
        case 3
            % 3D data maps indexed by dop (day of period)
            LAI  = LAI_map(:, :, dop);
            SAI  = SAI_map(:, :, dop);
            RD   = RD_map(:, :, dop);
            h_c  = h_c_map(:, :, dop);
            d    = d_map(:, :, dop);
            z_0m = z0M_map(:, :, dop);
            z_0H = z0H_map(:, :, dop);
            Alb  = Alb_map(:, :, dop);
    
            % Same photochemical accumulation flags with spatial maps
            PAC_flag = (TT_data(:, :, dop) > A_start) & (TT_data(:, :, dop) < A_end);
            PAC_flag_spring = (TT_data(:, :, dop) > A_start_spring) & (TT_data(:, :, dop) < A_end_spring);
            PAC_flag_medit  = (TT_data(:, :, dop) > A_start_medit) & (TT_data(:, :, dop) < A_end_medit);
    
            % Adjust PAC flags for regional climate classes
            PAC_flag(whclim_map == 1.5)  = PAC_flag_spring(whclim_map == 1.5);
            PAC_flag(whclim_map == -0.5) = PAC_flag_medit(whclim_map == -0.5);
    
            TT = TT_data(:, :, dop);
    end
    
    % Fix measurement height errors: if z_m_O3 is below canopy height, replace with average
    z_m_O3_tmp(z_m_O3_tmp <= h_c) = mean(z_m_O3_tmp, 'all');
    
    % Calculate maximum water holding capacity (mm) based on soil water content and root depth
    MWHC = Lm3 .* RD;
    
    % Calculate root surplus water due to root depth growth since previous timestep
    root_surplus = Lm3 .* (RD - RD_prev);
    
    % Update previous root depth for next timestep
    RD_prev = RD;
    
    % Display progress in console
    disp(['processing ' datestr(date_array(i_timestep))]);
    
    %% Calculated variables, conversions, and hourly variables -------------------
    
    % Convert temperature from Celsius to Kelvin
    tas_K = tas_tmp + config.const.T_0;
    
    % Calculate saturation vapor pressure (kPa) and its slope (delta)
    e_sat_kPa = T_M_e_sat(tas_tmp);
    delta     = T_M_delta(tas_tmp);
    
    % Calculate Vapor Pressure Deficit (VPD) in kPa, ensuring it’s non-negative
    VPD_kPa = max(0, e_sat_kPa - e_tmp);
    
    % Convert surface ozone from ppb to µg/m³ at measurement height
    sfo3_zm_ugm3 = ppb_2_ugm3(tas_tmp, config.const.M_O3, sfo3_tmp, ps_tmp);
    
    %% LIGHT Module ---------------------------------------------------------------
    
    % Compute photosynthetically active radiation (PAR) components and flags
    [PAR_tot_act, PAR_dir, PAR_diff, ...
        PAR_sun, PAR_shad, light_flag, ...
        PPFD_sun, PPFD_shad] = ...
        module_light_v1(B_tmp, rsds_tmp, ps_tmp, LAI, ...
        n_row, n_col, config);
    
    % Reset VPD sum when transitioning from dark to light conditions (sunrise)
    sumVPD(light_flag & ~light_flag_prev) = 0;
    light_flag_prev = light_flag;
    
    %% ATMOSPHERE Module - Part 1 ------------------------------------------------
    
    % Calculate sunlit/shaded LAI, soil shortwave radiation, radiation components,
    % and atmospheric properties needed for further modeling
    [LAI_sun, LAI_shad, N, Q_sw_soil, ...
        Rn_tmp, G_Wm2, Rn_soil, ...
        q_adim, c_p_JkgK, rho_air_kgm3, ...
        Sc_O3, Sc_W, Sc_H] = ...
        module_atmosphere_part_one(rsds_tmp, tas_tmp, tas_K, e_tmp, ...
        LAI, SAI, B_tmp, Alb, ps_tmp, Rn_tmp, ...
        config);
    
    %% SOIL Module ---------------------------------------------------------------
    
    % Compute soil water content, plant available water (PAW), available water capacity (AWC),
    % and update soil moisture states including root water surplus
    [W_cont, SWC_tmp, PAW_tmp, AWC_prev, AWC] = ...
        module_soil(AWC_prev, W_in_prev, evap_W_prev, tras_W_prev, root_surplus, ...
        W_in_prev, MWHC, config, par_g_stom_winter, WP, FC);
    
    %% STOMATAL Module -----------------------------------------------------------
    
    % Calculate stomatal conductance fractions for sun and shade leaves,
    % update VPD sums, stomatal response functions, and other physiological variables
    [g_s_frac_sun, g_s_frac_shad, sumVPD, jarvis_func, Uddling, ...
        g_s_W_sun_ms, g_s_O3_sun_ms, g_s_W_shad_ms, g_s_O3_shad_ms, ...
        g_s_act_W_ms, g_s_act_O3_ms, g_s_frac_sun_prev, g_s_frac_shad_prev] = ...
        module_stomatal_wheat_v1(tas_tmp, VPD_kPa, W_cont, PPFD_sun, PPFD_shad, co2_tmp, ...
        POD0_prev, ps_tmp, sumVPD, ...
        g_s_frac_sun_prev, g_s_frac_shad_prev, g_max_W_winter, g_max_O3_winter, g_max_W_medite, g_max_O3_medite, ...
        LAI_sun, LAI_shad, LAI, ...
        TT, ...
        par_g_stom_winter, par_g_stom_spring, par_g_stom_medite, whclim_map, config);
    
    %% ATMOSPHERE Module - Part 2 ------------------------------------------------
    
    % Complete atmospheric flux and concentration calculations, including turbulent transfer,
    % ozone and water vapor resistances, deposition velocities, and leaf-level variables
    [H_hour, u_star_hour, L_hour, u_hc, ...
        R_aH_hc_zmO3, R_bO3, R_surf_O3, R_aH_dz0m_zmT, R_bW, R_bH, R_inc, ...
        r_s_act_O3, r_s_act_W, rc, rbO3, g_ext, ...
        O3_h_c_ugm3, O3_h_c_ppb, O3_2m_ppb, O3_3m_ppb, F_O3, F_O3_ppb, F_totDDIM_O3_ppb, ...
        PODY, PODY_prev, POD0, POD0_prev, ...
        tleaf_C, VPD_leaf_kPa, e_sat_tleaf_kPa, ...
        tras_W, E_soil, E_wet, evap_W, evaptras_W, tras_W_prev, evap_W_prev, ...
        S_c, W_in, W_in_prev] = ...
        module_atmosphere_part_two(z_m_O3_tmp, tas_tmp, tas_K, e_sat_kPa, e_tmp, VPD_kPa, delta, ...
        sfcWind_tmp, pr_tmp, sfo3_zm_ugm3, H_tmp, ...
        light_flag, Rn_tmp, Rn_soil, G_Wm2, N, ...
        ps_tmp, rho_air_kgm3, c_p_JkgK, Sc_W, Sc_H, Sc_O3, ...
        d, z_0m, LAI, SAI, h_c, PAC_flag, ...
        g_s_act_O3_ms, g_s_act_W_ms, g_s_O3_sun_ms, PODY_prev, POD0_prev, ...
        AWC, S_c, W_in_prev, PAW_tmp, ...
        config, par_g_stom_winter);

    %% Accumulation of variables during the accumulation period (only light hours) ----------------
    
    % Increment count of timesteps within the accumulation period under light conditions
    n_timestep_PAClight = n_timestep_PAClight + double(PAC_flag & light_flag);
    
    % Cumulate Jarvis function components for temperature, soil moisture, and VPD
    jarvis_out_temp(PAC_flag & light_flag) = jarvis_func.ftemp(PAC_flag & light_flag) + jarvis_out_temp(PAC_flag & light_flag);
    jarvis_out_soil(PAC_flag & light_flag) = jarvis_func.fsoil(PAC_flag & light_flag) + jarvis_out_soil(PAC_flag & light_flag);
    jarvis_out_VPD(PAC_flag & light_flag)  = jarvis_func.fVPD (PAC_flag & light_flag) + jarvis_out_VPD(PAC_flag & light_flag);
    
    % Accumulate surface ozone concentration (ppb), temperature, VPD, and plant available water
    sfo3_ppb_out(PAC_flag & light_flag) = sfo3_ppb_out(PAC_flag & light_flag) + sfo3_tmp(PAC_flag & light_flag);
    T_out(PAC_flag & light_flag)         = T_out(PAC_flag & light_flag) + tas_tmp(PAC_flag & light_flag);
    VPD_out(PAC_flag & light_flag)       = VPD_out(PAC_flag & light_flag) + VPD_kPa(PAC_flag & light_flag);
    PAW_out(PAC_flag & light_flag)       = PAW_out(PAC_flag & light_flag) + PAW_tmp(PAC_flag & light_flag);
    
    % Save 2m and 3m ozone concentrations if running historical experiment
    if strcmp(config.experiment,'historical')
        o3_2m_output_ppb(:,:,i_timestep) = single(O3_2m_ppb);
        o3_3m_output_ppb(:,:,i_timestep) = single(O3_3m_ppb);
    end
    
    % Accumulate ozone difference (surface minus canopy height)
    sfo3_diff_ppb(PAC_flag & light_flag) = nansum(cat(2, ...
        sfo3_tmp(PAC_flag & light_flag) - O3_h_c_ppb(PAC_flag & light_flag), ...
        sfo3_diff_ppb(PAC_flag & light_flag)), 2);
    
    % Track the earliest day of accumulation period for flagged cells
    start_day_PAC(PAC_flag & light_flag) = min(start_day_PAC(PAC_flag & light_flag), dop);
    
    %% Save POD maps each year when maturity is reached -----------------------------------
    
    % Identify locations where maturity occurs: LAI drops from >0 to 0
    flag_maturity = (LAI_prev > 0) & (LAI == 0);
    
    % Save final values at maturity for POD and Jarvis functions (averaged over accumulation)
    PODY_final(flag_maturity)         = PODY(flag_maturity);
    POD0_final(flag_maturity)         = POD0(flag_maturity);
    jarvis_out_temp_final(flag_maturity) = jarvis_out_temp(flag_maturity) ./ n_timestep_PAClight(flag_maturity);
    jarvis_out_soil_final(flag_maturity) = jarvis_out_soil(flag_maturity) ./ n_timestep_PAClight(flag_maturity);
    jarvis_out_VPD_final(flag_maturity)  = jarvis_out_VPD(flag_maturity) ./ n_timestep_PAClight(flag_maturity);
    sfo3_ppb_out_final(flag_maturity)     = sfo3_ppb_out(flag_maturity) ./ n_timestep_PAClight(flag_maturity);
    sfo3_diff_ppb_final(flag_maturity)    = sfo3_diff_ppb(flag_maturity) ./ n_timestep_PAClight(flag_maturity);
    T_out_final(flag_maturity)             = T_out(flag_maturity) ./ n_timestep_PAClight(flag_maturity);
    VPD_out_final(flag_maturity)           = VPD_out(flag_maturity) ./ n_timestep_PAClight(flag_maturity);
    PAW_out_final(flag_maturity)           = PAW_out(flag_maturity) ./ n_timestep_PAClight(flag_maturity);
    n_day_PAC_final(flag_maturity)         = n_timestep_PAClight(flag_maturity) / (24 / config.timestep);
    start_day_PAC_final(flag_maturity)     = start_day_PAC(flag_maturity);
    
    % Reset variables after maturity to prepare for next accumulation period
    POD0_prev(flag_maturity)         = 0;
    PODY_prev(flag_maturity)         = 0;
    n_timestep_PAClight(flag_maturity) = 0;
    jarvis_out_temp(flag_maturity)  = 0;
    jarvis_out_soil(flag_maturity)  = 0;
    jarvis_out_VPD(flag_maturity)   = 0;
    PAW_out(flag_maturity)          = 0;
    sfo3_ppb_out(flag_maturity)     = 0;
    sfo3_diff_ppb(flag_maturity)    = 0;
    start_day_PAC(flag_maturity)    = nan;
    
    % Update previous LAI for next timestep comparison
    LAI_prev = LAI;
    
    %% Save yearly POD results to file ------------------------------------------------------
    
    if (~strcmp(year(date_array(min(i_timestep+1,length(date_array)))), year(date_array(i_timestep)))) || ...
            i_timestep == length(date_array))
        
        filename = fullfile(folder_save, ...
            ['POD' int2str(config.plant.Y) '_' int2str(year(date_array(i_timestep))) '.mat']);
        
        save(filename, 'POD0_final', 'PODY_final', 'jarvis_out_temp_final', ...
            'jarvis_out_soil_final', 'jarvis_out_VPD_final', 'sfo3_ppb_out_final', ...
            'sfo3_diff_ppb_final', 'n_day_PAC_final', 'start_day_PAC_final', ...
            'PAW_out_final', 'T_out_final', 'VPD_out_final');
    end
    
    %% Diagnostic table output --------------------------------------------------------
    
    if config.table_out_option
        new_row = cell(1, numel(var_names_array));
        
        for i_var = 1:numel(var_names_array)
            new_row{1, i_var} = eval(var_names_array{i_var});
        end
        
        diag_table(i_timestep, :) = new_row;
    end
end

% Display elapsed time for deposition model in minutes
disp(['deposition model: ' num2str(toc/60) ' min'])

%% Save initial settings for record keeping
filename = fullfile(folder_save, 'ini_settings.mat');
save(filename, 'config');

%% For the 'historical' experiment: save recursive state variables
% These variables are needed for initializing other experiments and to
% validate o3 concentrations at 3 and 2 meters during the period 2000-2014
% using the TOAR database

if strcmp('historical', config.experiment)
    % Save previous year’s accumulation and Jarvis function variables
    filename = fullfile(folder_main, folder_data, ...
        [model_name '_historical_prevData_' config.filter_hem '_' int2str(year(date_array(i_timestep))) '.mat']);
    save(filename, 'POD0_prev', 'PODY_prev', 'n_timestep_PAClight', ...
        'jarvis_out_temp', 'jarvis_out_soil', 'jarvis_out_VPD', ...
        'sfo3_ppb_out', 'sfo3_diff_ppb');
    
    % Save 3m ozone concentration data (large file, use v7.3 format)
    filename = fullfile(folder_save, 'O3_3m.mat');
    save(filename, 'o3_3m_output_ppb', 'lat', 'lon', 'time_tgt_days', 'date_array', '-v7.3');
    
    % Save 2m ozone concentration data (large file, use v7.3 format)
    filename = fullfile(folder_save, 'O3_2m.mat');
    save(filename, 'o3_2m_output_ppb', 'lat', 'lon', 'time_tgt_days', 'date_array', '-v7.3');
end

%% Optionally save diagnostic data to Excel file
if config.table_out_option
    filename = fullfile(folder_save, ['diagnostic_data' tmp_suffix '.xlsx']); % Construct filename with suffix
    writetable(diag_table, filename, 'Sheet', 1);
end
