function Rsurf = R_surf(LAI, SAI, r_stom, r_cut_O3, R_inc, R_soil)
%R_SURF Calculates the total surface resistance.
%
%   INPUTS:
%       LAI       - Leaf Area Index (dimensionless)
%       SAI       - Surface Area Index (dimensionless)
%       r_stom    - Stomatal resistance, scaled by LAI and considering sunlit and shaded fractions [s/m]
%       r_cut_O3  - Cuticular resistance including ozone effects [s/m]
%       R_inc     - In-canopy aerodynamic resistance [s/m]
%       R_soil    - Soil surface resistance [s/m]
%
%   OUTPUT:
%       Rsurf     - Overall surface resistance [s/m]
%
% Author: PR Guaita, 2025

Rsurf = (LAI ./ r_stom + SAI ./ r_cut_O3 + 1 ./ (R_inc + R_soil)) .^ (-1);

end
