function e_sat = T_M_e_sat(tas)
%T_M_E_SAT Calculates saturation vapor pressure using the Tetens-Murray equation.
%
%   This function computes the saturation vapor pressure (in kPa) based on air
%   temperature using the Tetens-Murray formula, valid for both positive and
%   negative temperatures.
%
%   INPUT:
%       tas - Air temperature (°C), scalar or matrix
%
%   OUTPUT:
%       e_sat - Saturation vapor pressure (kPa), same size as tas
%
%   References:
%       - Tetens (1930)
%       - Murray (1967)
%
% Author: PR Guaita, 2025

e_sat = zeros(size(tas));

% For temperatures above 0 °C
e_sat(tas > 0) = 0.611 * exp((17.269 * tas(tas > 0)) ./ (tas(tas > 0) + 237.3));

% For temperatures at or below 0 °C
e_sat(tas <= 0) = 0.611 * exp((21.875 * tas(tas <= 0)) ./ (tas(tas <= 0) + 265.5));

end
