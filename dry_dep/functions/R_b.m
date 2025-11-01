function R_b_Chim = R_b(config_const, u_star, Sc)
%R_B Calculates the bulk resistance of the leaf boundary laminar sub-layer.
%
%   Reference: Hick et al. (1987)
%
%   INPUTS:
%       config_const - Struct containing constants (e.g., von Karman constant and Prandtl number)
%       u_star       - Friction velocity [m/s]
%       Sc           - Schmidt number (dimensionless; depends on the molecule considered)
%
%   OUTPUT:
%       R_b_Chim     - Bulk resistance [s/m]
%
% Author: PR Guaita, 2025

kar = config_const.kar;
Pr  = config_const.Pr;

R_b_Chim = (2 ./ (kar .* u_star)) .* (Sc / Pr).^(2/3);

end
