function Tras_H2O = tras_stom_mms(T_C, rho, c_p, ...
    VPD_leaf_kPa, ...
    R_aH_dz0m_zmT, R_bH, ...
    delta, Rn, G, ...
    n, gamma, R_bW, ...
    r_s_act_W, LAI)
%TRAS_STOM_MMS Computes transpiration through stomata [mm s^-1]
% using the Penman-Monteith equation reformulated by Gerosa (2007, 2011)
% for stomatal conductance modelling.
%
%   INPUT:
%     T_C              - Air temperature [°C]
%     rho              - Air density [kg m^-3]
%     c_p              - Specific heat capacity of air at constant pressure [J kg^-1 K^-1]
%     VPD_leaf_kPa     - Vapor pressure deficit between leaf and air [kPa]
%     R_aH_dz0m_zmT    - Aerodynamic resistance between d+z0 and measurement height [s m^-1]
%     R_bH             - Quasi-laminar boundary layer resistance for heat [s m^-1]
%     delta            - Slope of saturation vapor pressure curve [kPa K^-1]
%     Rn               - Net radiation [W m^-2]
%     G                - Soil heat flux [W m^-2]
%     n                - 1 if vegetation is amphistomatous; 2 if hypostomatous
%     gamma            - Psychrometric constant [kPa K^-1]
%     R_bW             - Quasi-laminar boundary layer resistance for water [s m^-1]
%     r_s_act_W        - Actual stomatal resistance to water transfer [s m^-1]
%     LAI              - Leaf area index [m^2 m^-2]
%
%   OUTPUT:
%     Tras_H2O         - Instantaneous transpiration through stomata [mm s^-1]
%
%   Note:
%     r_s_act_W is already scaled by the fraction of sunlit and shaded leaves.
%     If VPD_leaf_kPa is negative or zero (i.e., e_sat < e), transpiration is forced to zero.
%
% Author: PR Guaita, 2025

% Latent heat of vaporization [J kg^-1], temperature dependent
lambda_J_kg = (-2733.33333333 * T_C + 2501333.33333);

% Compute transpiration flux [W m^-2] based on Penman-Monteith
Tras_H2O_Wm2 = ...
    (rho .* c_p .* (1000 * max(0, VPD_leaf_kPa)) ./ (R_aH_dz0m_zmT + R_bH) + ...
     1000 * delta .* (Rn - G)) ./ ...
    (1000 * delta + n * gamma .* (R_aH_dz0m_zmT + R_bW + r_s_act_W ./ LAI) ./ ...
    (R_aH_dz0m_zmT + R_bH));

% Convert from [W m^-2] to [mm s^-1]
Tras_H2O = Tras_H2O_Wm2 ./ lambda_J_kg;

% Safety check: warn if instantaneous transpiration exceeds 1 mm per hour
if any(Tras_H2O * 3600 > 1, 'all')
    disp('Warning! Transpiration exceeds 1 mm hr^-1')
end

end
