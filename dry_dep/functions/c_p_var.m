function c_p = c_p_var(c_p_std, q_adim)
%C_P_VAR Computes the specific heat capacity of moist air at constant pressure.
%
%   Calculates the effective specific heat capacity of air, accounting for 
%   the presence of water vapor.
%
%   INPUT:
%     c_p_std - Standard specific heat of dry air [J kg^-1 K^-1]
%     q_adim  - Specific humidity [dimensionless, kg H2O per kg moist air]
%
%   OUTPUT:
%     c_p     - Specific heat capacity of moist air [J kg^-1 K^-1]
%
%   Formula:
%     c_p = c_p_std * (1 + 0.84 * q_adim)
%
%   Note:
%     - The factor 0.84 accounts for the difference between the specific heat 
%       of water vapor (~1860 J kg^-1 K^-1) and dry air (~1005 J kg^-1 K^-1).
%
% Author: PR Guaita, 2025

c_p = c_p_std * (1 + 0.84 * q_adim);

end
