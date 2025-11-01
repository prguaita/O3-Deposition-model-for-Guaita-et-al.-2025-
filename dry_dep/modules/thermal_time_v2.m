function [TT_data] = thermal_time_v2(tas_C, tas_time, config, sow_date_map, TT_map_init)
%THERMAL_TIME_V2 Calculates accumulated thermal time starting from the sowing date each year.
%
% INPUTS:
%   tas_C         - Temperature data [°C], can be 2D (time) or 3D (lon x lat x time)
%   tas_time      - Time vector corresponding to tas_C data [days], fractional days allowed
%   config        - Configuration struct containing calendar options, starting year, etc.
%   sow_date_map  - Static map (lon x lat x years) with sowing day of year (DOY) for each location/year
%   TT_map_init   - Initial thermal time map (lon x lat), usually zeros or previous accumulation
%
% OUTPUT:
%   TT_data       - Thermal time accumulation map over time [lon x lat x days] or [days] if 1D data
%
% NOTES:
%   - Thermal time accumulates daily from sowing date resetting yearly
%   - Supports leap and fixed calendar years
%   - Supports 2D or 3D temperature inputs
%
% Author: PR Guaita, Aug 2023 (rewritten for clarity and documentation)

% Round down fractional days to obtain integer day indices for time series
time_serie_day = floor(tas_time);

% Calculate year and days per year for each timestep based on calendar option
switch config.input_model.opt_calendar
    case 'leap'
        % Compute calendar year for each timestep (supports leap years)
        time_serie_year = year(datetime(config.yr_start-1,12,31) + days(time_serie_day));
        days_in_year = 365 * ones(size(time_serie_year)) + leapyear(time_serie_year); % 366 for leap years
    case 'fix'
        % Fixed-length year calendar
        time_serie_year = config.yr_start + floor((time_serie_day - 1) ./ config.input_model.n_day_year);
        days_in_year = config.input_model.n_day_year * ones(size(time_serie_year));
end

% Pre-allocate thermal time data array: (lon x lat x unique days) or (unique days) for 2D data
TT_data = squeeze(nan([size(sow_date_map,1), size(sow_date_map,2), length(unique(time_serie_day))]));

% Flags to track sowing day resets for each grid cell at previous timestep
flag_sow_day_map_prev = false(size(sow_date_map,1), size(sow_date_map,2));

% Initialize cumulative thermal time map with initial values
TT_prev = TT_map_init;

% Temporary storage for accumulating temperature within a day
tmp_tas = zeros(size(sow_date_map,1), size(sow_date_map,2));

% Initialize arrays to flag new days and new years at each timestep
flag_new_day = zeros(size(tas_time));
flag_new_year = zeros(size(tas_time));

% Compute day and year change flags (1 if new day/year, else 0)
for i_timestep = 2:length(tas_time)
    flag_new_day(i_timestep) = time_serie_day(i_timestep) - time_serie_day(i_timestep-1);
    flag_new_year(i_timestep) = time_serie_year(i_timestep) - time_serie_year(i_timestep-1);
end

% Start with first year’s sowing date map slice
k_year = 1;
sow_date_map_year = sow_date_map(:,:,k_year);

% Loop over each timestep
for i_timestep = 1:length(tas_time)
    
    % If a new calendar year started, update the sowing date map for that year
    if flag_new_year(i_timestep)
        k_year = min(k_year + 1, size(sow_date_map, 3)); % avoid exceeding available years
        sow_date_map_year = sow_date_map(:,:,k_year);
    end
    
    % If a new day started
    if flag_new_day(i_timestep)
        % Reset temporary daily temperature accumulator
        tmp_tas = zeros(size(sow_date_map,1), size(sow_date_map,2));
        
        % Retrieve previous thermal time accumulation for the previous day
        switch ndims(tas_C)
            case 2
                TT_prev = TT_data(time_serie_day(i_timestep-1));
            case 3
                TT_prev = TT_data(:,:,time_serie_day(i_timestep-1));
        end
        
        % Identify locations where current day is sowing day (normalized by days_in_year)
        % Using modulo to handle day-of-year cycles, reset thermal time at sowing day
        flag_sow_day_map = (sow_date_map_year - 1 ./ days_in_year(i_timestep)) < ...
            mod(tas_time(i_timestep) ./ days_in_year(i_timestep), 1);
        
        % Reset thermal time at sowing day but only if it wasn't reset in the previous timestep
        TT_prev(flag_sow_day_map & ~flag_sow_day_map_prev) = 0;
        
        % Update previous sowing day flags
        flag_sow_day_map_prev = flag_sow_day_map;
    end
    
    % Accumulate thermal time within the day by adding current temperature weighted by timestep fraction
    switch ndims(tas_C)
        case 2
            % Compute weighting based on number of timesteps per day
            steps_per_day = find(floor(tas_time) > 1, 1) - 1;
            tmp_tas = tmp_tas + tas_C(i_timestep) / steps_per_day;
            
            % Add accumulated temperature (clipped at zero) to previous thermal time and store
            TT_data(time_serie_day(i_timestep)) = single(round(TT_prev + max(tmp_tas, 0, 'includemissing'), 1));
            
        case 3
            steps_per_day = find(floor(tas_time) > 1, 1) - 1;
            tmp_tas = tmp_tas + tas_C(:,:,i_timestep) / steps_per_day;
            
            TT_data(:,:,time_serie_day(i_timestep)) = single(round(TT_prev + max(tmp_tas, 0, 'includemissing'), 1));
    end
end

end
