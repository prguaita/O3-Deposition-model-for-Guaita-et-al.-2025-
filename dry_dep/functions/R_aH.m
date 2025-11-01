function R_a_zlow_zup = R_aH(z_low, z_up, d, kar, u_star, L)
%R_AH Calculates aerodynamic resistance between two heights using stability corrections.
%
%   INPUTS:
%       z_low   [m] Lower height in the calculation
%       z_up    [m] Upper height in the calculation
%       d       [m] Displacement height of the surface
%       kar     [--] Von Karman constant
%       u_star  [m/s] Friction velocity for the surface
%       L       [m] Monin-Obukhov length (stability parameter)
%
%   OUTPUT:
%       R_a_zlow_zup [s/m] Aerodynamic resistance between z_low and z_up
%
% Author: PR Guaita, 2025

R_a_zlow_zup = (1 ./ (kar * u_star)) .* ...
    (log((z_up - d) ./ (z_low - d)) - func_PSI_H((z_up - d) ./ L) + func_PSI_H((z_low - d) ./ L));

end
