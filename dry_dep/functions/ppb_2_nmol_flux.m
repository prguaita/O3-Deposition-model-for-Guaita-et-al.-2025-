function F_nmol = ppb_2_nmol_flux(F_ppb, R_univ, P, T_K)
%PPB_2_NMOL_FLUX Converts fluxes from ppb·m·s⁻¹ to nmol·m⁻²·s⁻¹
%
%   INPUTS:
%       F_ppb  - Flux in ppb·m·s⁻¹
%       R_univ - Universal gas constant [J mol⁻¹ K⁻¹]
%       P      - Pressure [Pa]
%       T_K    - Temperature [K]
%
%   OUTPUT:
%       F_nmol - Flux in nmol·m⁻²·s⁻¹
%
% Author: PR Guaita, 2025

F_nmol = F_ppb .* P ./ (R_univ * T_K);

end
