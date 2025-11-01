function [PAR_sun, PAR_sha] = PAR_leaf_fractions(PAR_dir, PAR_diff, LAI, B, gamma_leaf)
%PAR_LEAF_FRACTIONS Calculates PAR fractions on sunlit and shaded leaves.
%
%   This function estimates the photosynthetically active radiation (PAR) 
%   intercepted by sunlit and shaded leaves based on direct and diffuse PAR,
%   leaf area index (LAI), solar elevation angle, and leaf inclination angle.
%
% INPUTS:
%   PAR_dir     - Direct PAR component [µmol m⁻² s⁻¹]
%   PAR_diff    - Diffuse PAR component [µmol m⁻² s⁻¹]
%   LAI         - Leaf Area Index
%   B           - Solar elevation angle [degrees]
%   gamma_leaf  - Leaf inclination angle relative to the sun [degrees]
%
% OUTPUTS:
%   PAR_sun     - PAR on sunlit leaves [µmol m⁻² s⁻¹]
%   PAR_sha     - PAR on shaded leaves [µmol m⁻² s⁻¹]
%
% Notes:
%   - Based on Norman (1986)
%   - Leaf inclination is treated as a constant rather than a distribution.
%
% Author: PR Guaita, 2024

PAR_sha = PAR_diff .* exp(-0.5 .* LAI .^ 0.8) + 0.07 .* PAR_dir .* (1.1 - 0.1 .* LAI) .* exp(-sind(B));
PAR_sun = PAR_dir .* 0.8 .* cosd(gamma_leaf) ./ sind(B) + PAR_sha;

end
