function E_wet = evap_sur_wet_mms(T_C, rho, c_p, ...
    VPD_kPa, R_aH_dz0m_zmT, R_bH, ...
    delta, Rn, G, ...
    n, gamma, R_bW)
%EVAP_SUR_WET_MMS Computes evaporation from wet surfaces [mm s^-1]
% using the Penman-Monteith equation adapted for wet canopy or soil surfaces.
%
%   INPUT:
%     T_C              - Air temperature [°C]
%     rho              - Air density [kg m^-3]
%     c_p              - Specific heat capacity of air at constant pressure [J kg^-1 K^-1]
%     VPD_kPa          - Vapor pressure deficit at measurement height [kPa]
%     R_aH_dz0m_zmT    - Aerodynamic resistance between d+z0 and measurement height [s m^-1]
%     R_bH             - Quasi-laminar boundary layer resistance for heat [s m^-1]
%     delta            - Slope of the saturation vapor pressure curve [kPa K^-1]
%     Rn               - Net radiation at wet surface [W m^-2]
%     G                - Ground heat flux [W m^-2]
%     n                - Number of evaporating surfaces (1=amphistomatous, 2=hypostomatous)
%     gamma            - Psychrometric constant [kPa K^-1]
%     R_bW             - Quasi-laminar boundary layer resistance for water vapor [s m^-1]
%
%   OUTPUT:
%     E_wet            - Instantaneous evaporation from wet surfaces [mm s^-1]
%
%   Notes:
%     - Evaporation is forced to zero if calculation becomes negative.
%     - Latent heat of vaporization depends on temperature.
%
% Author: PR Guaita, 2025

% Latent heat of vaporization [J kg^-1], temperature dependent
lambda_J_kg = (-2733.33333333 * T_C + 2501333.33333);

% Compute latent heat flux [W m^-2]
LE_Wm2 = ...
    (rho .* c_p .* (VPD_kPa * 1000) ./ (R_aH_dz0m_zmT + R_bH) + ...
     1000 * delta .* (Rn - G)) ./ ...
    (1000 * delta + ...
     n * gamma .* (R_aH_dz0m_zmT + R_bW) ./ (R_aH_dz0m_zmT + R_bH));

% Convert from [W m^-2] to [mm s^-1] and force negative values to zero
E_wet = max(0, LE_Wm2 ./ lambda_J_kg);

end
