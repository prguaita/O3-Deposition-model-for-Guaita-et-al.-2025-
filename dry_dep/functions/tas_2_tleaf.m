function tleaf = tas_2_tleaf(tas, H, Ra_hc_zmT, rb_heat, rho, c_p)
%TAS_2_TLEAF Converts air temperature to leaf temperature.
%
%   This function estimates leaf temperature based on air temperature and
%   sensible heat flux, accounting for resistances between measurement heights.
%
%   INPUTS:
%       tas         - Air temperature at reference height (°C), typically 2 m
%       H           - Sensible heat flux (W m⁻²)
%       Ra_hc_zmT   - Aerodynamic resistance between air temperature measurement height 
%                     and canopy height (s m⁻¹)
%       rb_heat     - Boundary layer resistance of a single leaf (s m⁻¹)
%       rho         - Air density (kg m⁻³)
%       c_p         - Specific heat capacity of air at constant pressure (J kg⁻¹ K⁻¹)
%
%   OUTPUT:
%       tleaf       - Estimated leaf temperature (°C), capped at 6°C above air temperature
%
% Author: PR Guaita, 2025

tleaf = tas + H .* (Ra_hc_zmT + rb_heat) ./ (rho .* c_p);
tleaf = min(tleaf, tas + 6);  % Cap leaf temperature at max 6°C above air temperature

end
