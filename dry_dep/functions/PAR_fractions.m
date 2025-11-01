function [PAR_dir, PAR_diff] = PAR_fractions(PAR, P, B)
%PAR_FRACTIONS Splits total PAR into direct and diffuse components.
%
%   This function estimates the direct and diffuse Photosynthetically Active Radiation (PAR)
%   fractions based on total PAR, atmospheric pressure, and solar elevation angle.
%   The calculation follows Weiss and Norman (1985) methodology.
%
% INPUTS:
%   PAR - Total Photosynthetically Active Radiation [µmol m⁻² s⁻¹]
%   P   - Atmospheric pressure [kPa]
%   B   - Solar elevation angle [degrees]
%
% OUTPUTS:
%   PAR_dir  - Direct (beam) PAR component [µmol m⁻² s⁻¹]
%   PAR_diff - Diffuse PAR component [µmol m⁻² s⁻¹]
%
% Notes:
%   - Optical air mass is computed from the solar elevation angle.
%   - Potential direct and diffuse PAR are estimated using empirical relationships.
%   - Sky transmissivity adjusts the partitioning of PAR into direct and diffuse components.
%
% Author: PR Guaita, 2024

% Constants
P_0 = 101.325; % Standard atmospheric pressure at sea level [kPa]

% Calculate optical air mass (relative air mass)
m = 1 ./ sind(B);

% Estimate potential direct PAR (pPAR_dir) based on pressure and solar elevation
pPAR_dir = 600 .* exp(-0.185 .* (P ./ P_0) .* m) .* sind(B);

% Estimate potential diffuse PAR (pPAR_diff), constrained to be non-negative
pPAR_diff = 0.4 .* max(0, 600 - pPAR_dir) .* sind(B);

% Total potential PAR
pPAR_tot = pPAR_dir + pPAR_diff;

% Calculate sky transmissivity (bounded between 0.21 and 0.9)
ST = min(0.9, max(0.21, PAR ./ max(PAR, pPAR_tot)));

% Calculate fraction of direct PAR (adjusted by sky transmissivity)
fPAR_dir = (pPAR_dir ./ pPAR_tot) .* (1 - ((0.9 - ST) ./ 0.7).^(2/3));

% Fraction of diffuse PAR is complementary to direct fraction
fPAR_diff = 1 - fPAR_dir;

% Calculate direct and diffuse PAR components
PAR_dir = fPAR_dir .* PAR;
PAR_diff = fPAR_diff .* PAR;

end
