function conc_z = O3_z(conc_zm, R_a_z_zm, R_tot)
%O3_Z Calculates ozone concentration at height z.
%
%   The height z (target height) is implicit in the aerodynamic resistance
%   between z and the measurement height.
%
%   INPUTS:
%       conc_zm   - Ozone concentration at measurement height (units as given)
%       R_a_z_zm  - Aerodynamic resistance between height z and measurement height [s/m]
%       R_tot     - Total resistance [s/m]
%
%   OUTPUT:
%       conc_z    - Ozone concentration at height z (same units as conc_zm)
%
% Author: PR Guaita, 2025

conc_z = conc_zm .* (1 - R_a_z_zm ./ R_tot);

end
