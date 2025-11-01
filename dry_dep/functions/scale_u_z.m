function u_z = scale_u_z(z, u_star, kar, d, z_0, L)
%SCALE_U_Z Scales wind speed to a target height above the surface.
%
%   INPUTS:
%       z       [m]   Target height above ground where wind speed is computed
%       u_star  [m/s] Friction velocity for the considered surface
%       kar     [--]  Von Karman constant
%       d       [m]   Displacement height for the considered surface
%       z_0     [m]   Roughness length for the considered surface
%       L       [m]   Monin-Obukhov length for the considered surface
%
%   OUTPUT:
%       u_z     [m/s] Wind speed at height z above the considered surface
%
%   Note: Inputs u_star, d, and z_0 should correspond consistently to the same surface level (target or reference).
%
% Author: PR Guaita, 2025

u_z = (u_star / kar) .* (log((z - d) ./ z_0) - func_PSI_M((z - d) ./ L) + func_PSI_M(z_0 ./ L));

end
