function rho = air_density(q, T, P_z, config_const)
%AIR_DENSITY Calculates the density of moist air
%
%   INPUTS:
%     q          - Specific humidity [kg H2O per kg moist air]
%     T          - Air temperature [°C]
%     P_z        - Atmospheric pressure [kPa]
%     config_const - Structure containing constants:
%                   M_dry: Molar mass of dry air [g/mol]
%                   Re: Universal gas constant [J mol^-1 K^-1]
%
%   OUTPUT:
%     rho        - Air density [kg/m^3]
%
%   FORMULA:
%     rho = (M_dry * P_z) / (Re * (T + 273.15) * (1 + 0.61 * q))
%
%   Notes:
%     - Temperature converted to Kelvin inside the function.
%     - The factor (1 + 0.61*q) accounts for the effect of water vapor on air density.
%
% Author: PR Guaita, 2025

M_dry = config_const.M_dry;  % Molar mass of dry air [g/mol]
Re    = config_const.Re;     % Universal gas constant [J mol^-1 K^-1]

rho = (M_dry * P_z) ./ (Re * (T + 273.15) .* (1 + 0.61 * q));

end
