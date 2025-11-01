function delta = T_M_delta(T)
%T_M_DELTA Computes the slope of the saturation vapor pressure curve (delta)
% using the Tetens-Murray formula, based on temperature in °C.
%
%   INPUT:
%       T - Air temperature (°C), scalar
%
%   OUTPUT:
%       delta - Slope of the saturation vapor pressure curve (kPa °C⁻¹)
%
%   Note:
%       Uses Tetens-Murray coefficients separately for T >= 0 °C and T < 0 °C.
%
%   Verified by GG
%
% Author: PR Guaita, 2025

if T >= 0
    % For temperatures >= 0 °C
    e_exp = exp((17.27 * T) ./ (T + 237.3));
    delta = 0.61078 * e_exp .* (17.269 * 237.3) ./ ((T + 237.3).^2);
else
    % For temperatures < 0 °C
    e_exp = exp((21.875 * T) ./ (T + 265.5));
    delta = 0.61078 * e_exp .* (21.875 * 265.5) ./ ((T + 265.5).^2);
end

end
function delta = T_M_delta(T)
%T_M_DELTA Computes the slope of the saturation vapor pressure curve (delta)
% using the Tetens-Murray formula, based on temperature in °C.
%
%   INPUT:
%       T - Air temperature (°C), scalar
%
%   OUTPUT:
%       delta - Slope of the saturation vapor pressure curve (kPa °C-1)
%
%   Note:
%       Uses Tetens-Murray coefficients separately for T >= 0 °C and T < 0 °C.
%
%   Verified by GG
%
% Author: PR Guaita, 2025

if T >= 0
    % For temperatures >= 0 °C
    e_exp = exp((17.27 * T) ./ (T + 237.3));
    delta = 0.61078 * e_exp .* (17.269 * 237.3) ./ ((T + 237.3).^2);
else
    % For temperatures < 0 °C
    e_exp = exp((21.875 * T) ./ (T + 265.5));
    delta = 0.61078 * e_exp .* (21.875 * 265.5) ./ ((T + 265.5).^2);
end

end
