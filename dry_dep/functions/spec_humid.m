function q = spec_humid(e, P_z)
%SPEC_HUMID Computes the specific humidity.
%
%   Calculates the specific humidity of air given the vapor pressure and 
%   atmospheric pressure. The formula assumes pressures are in kPa.
%
%   INPUT:
%     e   - Water vapor pressure [kPa]
%     P_z - Atmospheric pressure at canopy or measurement height [kPa]
%
%   OUTPUT:
%     q   - Specific humidity [kg H2O per kg moist air]
%
%   Formula:
%     q = 0.622 * e / (P_z - e + 0.622 * e)
%
%   Note:
%     - When e and P_z are given in kPa, the output q is in [kg H2O / kg air].
%     - 0.622 is the ratio of molecular weight of water vapor (18.016 g/mol)
%       to dry air (28.964 g/mol).
%
% Author: PR Guaita, 2025

q = 0.622 * e ./ ((P_z - e) + 0.622 * e);

end
