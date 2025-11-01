function [LAI_sun, LAI_shad, N, Q_sw_soil, ...
          Rn_Wm2, G_Wm2, Rn_soil, ...
          q_adim, c_p_JkgK, rho_air_kgm3, ...
          Sc_O3, Sc_W, Sc_H] = ...
    module_atmosphere_part_one(Q_sw, T_C, T_K, e_kPa, ...
                               LAI, SAI, B, Alb, P_z_hour, Rn_Wm2, ...
                               config)
% MODULE_ATMOSPHERE_PART_ONE Computes atmospheric and radiation parameters
%
%   This module calculates sunlit and shaded Leaf Area Index (LAI),
%   cloud cover fraction, net radiation components, soil shortwave radiation,
%   specific humidity, specific heat capacity of air, air density, and Schmidt numbers.
%
% Inputs:
%   Q_sw        - Incoming shortwave radiation [W m^-2]
%   T_C         - Air temperature [°C]
%   T_K         - Air temperature [K]
%   e_kPa       - Vapor pressure [kPa]
%   LAI         - Leaf Area Index (dimensionless)
%   SAI         - Surface Area Index (dimensionless)
%   B           - Solar elevation (°)
%   Alb         - Surface albedo (dimensionless)
%   P_z_hour    - Atmospheric pressure at height z [kPa]
%   Rn_Wm2      - Net radiation at the canopy [W m^-2]; NaN if to be calculated
%   config      - Configuration struct with constants and parameters
%
% Outputs:
%   LAI_sun     - Sunlit Leaf Area Index
%   LAI_shad    - Shaded Leaf Area Index
%   N           - Cloud cover fraction (dimensionless)
%   Q_sw_soil   - Soil-absorbed shortwave radiation [W m^-2]
%   Rn_Wm2      - Net radiation at canopy [W m^-2]
%   G_Wm2       - Ground heat flux [W m^-2]
%   Rn_soil     - Net radiation at soil surface [W m^-2]
%   q_adim      - Specific humidity (dimensionless fraction)
%   c_p_JkgK    - Specific heat capacity of air at constant pressure [J kg^-1 K^-1]
%   rho_air_kgm3- Air density [kg m^-3]
%   Sc_O3       - Schmidt number for O3
%   Sc_W        - Schmidt number for water vapor (H2O)
%   Sc_H        - Schmidt number for heat

% Calculate sunlit and shaded leaf area index
[LAI_sun, LAI_shad] = LAI_ls(LAI, B);

% Estimate cloud cover fraction from shortwave radiation and beam fraction
N = CloudCover(Q_sw, B);

% Calculate net radiation if not provided
if isnan(Rn_Wm2)
    Rn_Wm2 = net_rad(T_K, config.const.alpha, ...
                       N, Alb, Q_sw);
end

% Ground heat flux approximated as 10% of net radiation
G_Wm2 = 0.1 * Rn_Wm2;

% Calculate soil absorbed shortwave radiation accounting for stem area index and canopy extinction coefficient
Q_sw_soil = soil_sw(Q_sw, config.const.Ka, SAI);

% Calculate net radiation at the soil surface
Rn_soil = net_rad(T_K, config.const.alpha, ...
                    N, Alb, Q_sw_soil);

% Calculate specific humidity (dimensionless)
q_adim = spec_humid(e_kPa, P_z_hour);

% Calculate specific heat capacity of air (J/kg/K) adjusted for humidity
c_p_JkgK = c_p_var(config.const.c_p_std, q_adim);

% Calculate air density [kg/m^3]
rho_air_kgm3 = air_density(q_adim, T_C, P_z_hour, config.const);

% Calculate Schmidt numbers for O3, water vapor, and heat
[Sc_O3, Sc_W, Sc_H] = Massman_Schmidt(T_C, P_z_hour, config.const);

end
