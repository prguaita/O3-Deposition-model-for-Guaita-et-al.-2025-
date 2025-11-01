function [g_s_frac_sun, g_s_frac_shad, jarvis_func] = ...
    jarvis_wheat_v1(par_g_stom, config, TT_data, POD0_prev, PPFD_sun, PPFD_shad, ...
    tas, VPD, soil, CO2_ppm)
%JARVIS_WHEAT_V1 Calculates relative stomatal conductance for wheat using the Jarvis model.
%
%   Computes stomatal conductance as a fraction of maximum conductance (0 to 1)
%   using a multiplicative Jarvis-type algorithm parameterized for winter wheat.
%   The function returns conductance fractions for sunlit and shaded canopy parts,
%   along with the intermediate multiplicative factors.
%
% INPUTS:
%   par_g_stom  - Struct with stomatal conductance parameters (e.g. T thresholds, light, VPD, soil)
%   config      - Struct with additional model configuration parameters
%   TT_data     - Thermal time (°C day) at current timestep (2D matrix)
%   POD0_prev   - Phytotoxic ozone dose above threshold (mmol m^-2) from previous timestep (2D matrix)
%   PPFD_sun    - Photosynthetic photon flux density in sunlit canopy [µmol m^-2 s^-1] (2D matrix)
%   PPFD_shad   - Photosynthetic photon flux density in shaded canopy [µmol m^-2 s^-1] (2D matrix)
%   tas         - Air temperature at surface [°C] (2D matrix)
%   VPD         - Vapor pressure deficit [kPa] (2D matrix)
%   soil        - Soil water content or plant available water fraction (0-1) (2D matrix)
%   CO2_ppm     - Atmospheric CO2 concentration [ppm] (scalar or 2D matrix)
%
% OUTPUTS:
%   g_s_frac_sun  - Fractional stomatal conductance for sunlit canopy [0-1] (2D matrix)
%   g_s_frac_shad - Fractional stomatal conductance for shaded canopy [0-1] (2D matrix)
%   jarvis_func   - Struct with intermediate Jarvis function components (fphen, fO3, flight_sun, flight_shad, ftemp, fVPD, fsoil)
%
% Author: PR Guaita - 2024

%% Get matrix dimensions
[n_row, n_col] = size(POD0_prev);

%% Extract parameters from structs
A_start = config.plant_TT.TT_A_start;
A_end   = config.plant_TT.TT_A_end;
p_a     = par_g_stom.p_a;
p_e     = par_g_stom.p_e;
p_one   = par_g_stom.p_one;
p_two   = par_g_stom.p_two;
p_three = par_g_stom.p_three;
f_min   = par_g_stom.f_min;

%% Phenology factor (fphen) based on thermal time

fphen = zeros(n_row, n_col);

case_one   = TT_data < A_start;
case_two   = (TT_data >= A_start) & (TT_data < p_one);
case_three = (TT_data >= p_one) & (TT_data < p_two);
case_four  = (TT_data >= p_two) & (TT_data < p_three);
case_five  = (TT_data >= p_three) & (TT_data < A_end);
case_six   = TT_data >= A_end;

fphen(case_one)   = TT_data(case_one) * p_a / A_start;
fphen(case_two)   = p_a + (1 - p_a) / (p_one - A_start) .* (TT_data(case_two) - A_start);
fphen(case_three) = 1;
fphen(case_four)  = 1 + (p_e - 1) / (p_three - p_two) .* (TT_data(case_four) - p_two);
fphen(case_five)  = p_e + (-p_e) / (A_end - p_three) .* (TT_data(case_five) - p_three);
fphen(case_six)   = 0;

%% Ozone damage factor (fO3)

if par_g_stom.f_ozone_option
    fO3 = (1 + (POD0_prev / 14).^8).^(-1);
else
    fO3 = ones(n_row, n_col);
end

%% Light response factors (flight) for sunlit and shaded canopy

flight_sun = 1 - exp(-par_g_stom.light_a * PPFD_sun);
flight_shad = 1 - exp(-par_g_stom.light_a * PPFD_shad);

%% Temperature response factor (ftemp)

T_max = par_g_stom.T_max;
T_opt = par_g_stom.T_opt;
T_min = par_g_stom.T_min;
bt = (T_max - T_opt) / (T_opt - T_min);

ftemp = nan(n_row, n_col);
% Set minimum conductance outside temperature range
ftemp(tas < T_min | tas > T_max) = f_min;

% Compute ftemp within optimal temperature range
valid_temp = (tas >= T_min) & (tas <= T_max);
ftemp(valid_temp) = max(f_min, ...
    ((tas(valid_temp) - T_min) / (T_opt - T_min)) .* ...
    (((T_max - tas(valid_temp)) / (T_max - T_opt)).^bt));

%% Vapor pressure deficit response factor (fVPD)

VPD_min = par_g_stom.VPD_min * ones(n_row, n_col);
VPD_max = par_g_stom.VPD_max * ones(n_row, n_col);

fVPD = min(1, max(f_min, ...
    (1 - f_min) .* (VPD_min - VPD) ./ (VPD_min - VPD_max) + f_min));

%% Soil moisture response factor (fsoil)

fsoil = zeros(n_row, n_col);

switch par_g_stom.soil_option
    case 'PAW'
        PAW_t = 0.01 * par_g_stom.PAW_t;
        fsoil(soil < PAW_t) = 1 + (soil(soil < PAW_t) - PAW_t) / PAW_t;
        fsoil(soil >= PAW_t) = 1;

    case 'SWC'
        SWC_min = 0.01 * par_g_stom.SWC_min;
        SWC_max = 0.01 * par_g_stom.SWC_max;

        in_range = (soil >= SWC_min) & (soil < SWC_max);
        fsoil(in_range) = max(f_min, min(1, ...
            f_min + (1 - f_min) .* (soil(in_range) - SWC_min) ./ (SWC_max - SWC_min)));

        fsoil(soil >= SWC_max) = 1;
        fsoil(soil < SWC_min) = f_min;

    case 'FC'
        fsoil = ones(n_row, n_col);

    otherwise
        error('Invalid "soil_option" specified in par_g_stom.');
end

%% CO2 response factor (fCO2)

fmin_CO2 = par_g_stom.fmin_CO2;
a_CO2 = par_g_stom.a_CO2;
CO2_ref = par_g_stom.CO2_ref;

if config.dep.CO2_effect
    fCO2 = fmin_CO2 + (1 - fmin_CO2) .* exp(-a_CO2 * (CO2_ppm - CO2_ref) / CO2_ref);
else
    fCO2 = 1;
end

%% Calculate final stomatal conductance fractions for sunlit and shaded leaves

g_s_frac_sun = min(fphen, fO3) .* flight_sun .* max(f_min, (ftemp .* fVPD .* fsoil)) * fCO2;
g_s_frac_shad = min(fphen, fO3) .* flight_shad .* max(f_min, (ftemp .* fVPD .* fsoil)) * fCO2;

%% Return intermediate Jarvis function components

jarvis_func = struct(...
    'fphen', fphen, ...
    'fO3', fO3, ...
    'flight_sun', flight_sun, ...
    'flight_shad', flight_shad, ...
    'ftemp', ftemp, ...
    'fVPD', fVPD, ...
    'fsoil', fsoil);

end
