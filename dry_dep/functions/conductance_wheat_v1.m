function [g_s_frac_sun, g_s_frac_shad, sumVPD, jarvis_func, Uddling] = ...
    conductance_wheat_v1(TT_data, par_g_stom_winter, par_g_stom_spring, par_g_stom_medite, ...
    whclim_map, config, POD0_prev, PPFD_sun, PPFD_shad, tas, VPD, PAW, CO2_ppm, ...
    sumVPD, g_s_frac_sun_prev, g_s_frac_shad_prev)
%CONDUCTANCE_WHEAT_V1 Calculates stomatal conductance fractions for wheat types.
%
%   This module computes the stomatal conductance as a fraction of maximum 
%   conductance (value between 0 and 1) for different wheat parametrizations 
%   (winter, spring, Mediterranean). The Uddling condition, which modifies 
%   conductance under high cumulative VPD stress, can optionally be applied.
%
% INPUTS:
%   TT_data            - Thermal time map (2D, lon x lat) at current timeframe
%   par_g_stom_winter  - Parameters struct for winter wheat stomatal conductance
%   par_g_stom_spring  - Parameters struct for spring wheat stomatal conductance
%   par_g_stom_medite  - Parameters struct for Mediterranean wheat stomatal conductance
%   whclim_map         - Wheat climate map (2D lon x lat), with codes indicating wheat type
%                        (1.5 = spring wheat, -0.5 or 0.5 = Mediterranean)
%   config             - Configuration struct with model parameters
%   POD0_prev          - Phytotoxic ozone dose above threshold at previous timestep [mmol m^-2]
%   PPFD_sun           - Photosynthetic photon flux density in sunlit canopy [µmol m^-2 s^-1]
%   PPFD_shad          - Photosynthetic photon flux density in shaded canopy [µmol m^-2 s^-1]
%   tas                - Air temperature at surface [°C]
%   VPD                - Vapor pressure deficit [kPa]
%   PAW                - Plant available water fraction (0 to 1)
%   CO2_ppm            - Atmospheric CO2 concentration [ppm]
%   sumVPD             - Cumulative VPD sum from sunrise [kPa]
%   g_s_frac_sun_prev  - Stomatal conductance fraction in sunlit leaves at previous timestep [0-1]
%   g_s_frac_shad_prev - Stomatal conductance fraction in shaded leaves at previous timestep [0-1]
%
% OUTPUTS:
%   g_s_frac_sun       - Updated stomatal conductance fraction for sunlit leaves [0-1]
%   g_s_frac_shad      - Updated stomatal conductance fraction for shaded leaves [0-1]
%   sumVPD             - Updated cumulative VPD sum from sunrise [kPa]
%   jarvis_func        - Struct containing intermediate Jarvis model factors (phenology, ozone, light, temperature, VPD, soil)
%   Uddling            - Logical map indicating grid cells where Uddling condition is active (true = conductance reduced)
%
% Author: PR Guaita - 2023

%% Identify wheat type flags based on climate map codes
flag_spring_wheat = (whclim_map == 1.5);
flag_medite_wheat = (whclim_map == -0.5) | (whclim_map == 0.5);

%% Calculate stomatal conductance fractions using Jarvis model for each wheat type

% Winter wheat parametrization (default)
[g_s_frac_sun, g_s_frac_shad, jarvis_func] = jarvis_wheat_v1(...
    par_g_stom_winter, config, TT_data, POD0_prev, PPFD_sun, PPFD_shad, tas, VPD, PAW, CO2_ppm);

% Spring wheat parametrization (temporary results)
[g_s_frac_sun_tmp, g_s_frac_shad_tmp, jarvis_func_tmp] = jarvis_wheat_v1(...
    par_g_stom_spring, config, TT_data, POD0_prev, PPFD_sun, PPFD_shad, tas, VPD, PAW, CO2_ppm);

% Replace winter wheat results with spring wheat where applicable
g_s_frac_sun(flag_spring_wheat) = g_s_frac_sun_tmp(flag_spring_wheat);
g_s_frac_shad(flag_spring_wheat) = g_s_frac_shad_tmp(flag_spring_wheat);

% Update Jarvis function factors for spring wheat locations
fields = fieldnames(jarvis_func);
for f = 1:numel(fields)
    field = fields{f};
    jarvis_func.(field)(flag_spring_wheat) = jarvis_func_tmp.(field)(flag_spring_wheat);
end

% Mediterranean wheat parametrization (temporary results)
[g_s_frac_sun_tmp, g_s_frac_shad_tmp, jarvis_func_tmp] = jarvis_wheat_v1(...
    par_g_stom_medite, config, TT_data, POD0_prev, PPFD_sun, PPFD_shad, tas, VPD, PAW, CO2_ppm);

% Replace results with Mediterranean wheat where applicable
g_s_frac_sun(flag_medite_wheat) = g_s_frac_sun_tmp(flag_medite_wheat);
g_s_frac_shad(flag_medite_wheat) = g_s_frac_shad_tmp(flag_medite_wheat);

% Update Jarvis function factors for Mediterranean wheat locations
for f = 1:numel(fields)
    field = fields{f};
    jarvis_func.(field)(flag_medite_wheat) = jarvis_func_tmp.(field)(flag_medite_wheat);
end

% Get output map size
[n_row, n_col] = size(g_s_frac_sun);

%% Apply Uddling condition to reduce conductance under high VPD stress
if par_g_stom_winter.uddling_option
    % Update cumulative VPD sum from sunrise
    sumVPD = sumVPD + VPD;
    
    % Initialize logical map for Uddling condition (conductance reduction)
    Uddling = false(n_row, n_col);
    
    % Mark grid cells where cumulative VPD exceeds critical threshold
    Uddling(sumVPD > par_g_stom_winter.sum_VPD_crit) = true;
    
    % Uddling reduction factor
    k_uddling = par_g_stom_winter.k_uddling;
    
    % Reduce stomatal conductance fractions proportionally in stressed cells
    g_s_frac_sun(Uddling) = min(g_s_frac_sun(Uddling), k_uddling * g_s_frac_sun_prev(Uddling));
    g_s_frac_shad(Uddling) = min(g_s_frac_shad(Uddling), k_uddling * g_s_frac_shad_prev(Uddling));
    
    % For cells not exceeding threshold, conductance remains unchanged
    % (This step is actually redundant but kept for clarity)
    g_s_frac_sun(~Uddling) = g_s_frac_sun(~Uddling);
    g_s_frac_shad(~Uddling) = g_s_frac_shad(~Uddling);
else
    % If Uddling condition inactive, initialize logical map to false
    Uddling = false(n_row, n_col);
end

end
