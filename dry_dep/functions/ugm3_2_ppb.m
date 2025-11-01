function conc_ppb = ugm3_2_ppb(T_K, molar_mass, conc_ugm3, P_z, R)
%UGM3_2_PPB Converts concentration from [µg/m³] to ppb
%
%   INPUTS:
%       T_K        - Temperature [K]
%       molar_mass - Molar mass of the gas [g/mol]
%       conc_ugm3  - Gas concentration [µg/m³]
%       P_z        - Atmospheric pressure at height z [kPa]
%       R          - Universal gas constant [J/(mol·K)]
%
%   OUTPUT:
%       conc_ppb   - Gas concentration in parts per billion (ppb)
%
% Author: PR Guaita, 2025

conc_ppb = conc_ugm3 .* (R * T_K) * 1000 ./ (molar_mass .* (P_z * 1000));

end
