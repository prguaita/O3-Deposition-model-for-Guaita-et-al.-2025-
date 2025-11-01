function E_soil = evap_soil_mms(T_C, rho, ...
    c_p, VPD_kPa, ...
    R_aH_dz0m_zmT, R_bH, R_inc, ...
    delta, Rn_soil, G, n, ...
    gamma, R_bW, R_soil)
%EVAP_SOIL_MMS Computes soil evaporation [mm s^-1]
% using the Penman-Monteith equation adapted for soil evaporation.
%
%   INPUT:
%     T_C              - Air temperature [°C]
%     rho              - Air density [kg m^-3]
%     c_p              - Specific heat capacity of air at constant pressure [J kg^-1 K^-1]
%     VPD_kPa          - Vapor pressure deficit at measurement height (not at soil surface) [kPa]
%     R_aH_dz0m_zmT    - Aerodynamic resistance between d+z0 and measurement height [s m^-1]
%     R_bH             - Quasi-laminar boundary layer resistance for heat [s m^-1]
%     R_inc            - Intra-canopy resistance [s m^-1]
%     delta            - Slope of saturation vapor pressure curve [kPa K^-1]
%     Rn_soil          - Net radiation at soil surface [W m^-2]
%     G                - Soil heat flux [W m^-2]
%     n                - 2 for soil (number of evaporating surfaces)
%     gamma            - Psychrometric constant [kPa K^-1]
%     R_bW             - Quasi-laminar boundary layer resistance for water [s m^-1]
%     R_soil           - Surface resistance of the soil [s m^-1]
%
%   OUTPUT:
%     E_soil           - Instantaneous soil evaporation [mm s^-1]
%
%   Notes:
%     - Evaporation is forced to zero when the calculation becomes negative
%       (e.g., at night due to negative VPD or energy imbalance).
%     - Latent heat of vaporization depends on temperature.
%
% Author: PR Guaita, 2025

% Latent heat of vaporization [J kg^-1], temperature dependent
lambda_J_kg = (-2733.33333333 * T_C + 2501333.33333);

% Compute latent heat flux [W m^-2] using Penman-Monteith
LE_Wm2 = ...
    (rho .* c_p .* (1000 * max(0, VPD_kPa)) ./ (R_aH_dz0m_zmT + R_bH + R_inc) + ...
     1000 * delta .* (Rn_soil - G)) ./ ...
    (1000 * delta + ...
     n * gamma .* (R_aH_dz0m_zmT + R_bW + R_inc + R_soil) ./ ...
    (R_aH_dz0m_zmT + R_bH + R_inc));

% Convert from [W m^-2] to [mm s^-1] and force negative values to zero
E_soil = max(0, LE_Wm2 ./ lambda_J_kg);

end
