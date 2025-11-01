function Fs_O3 = FsO3_MM(O3, g_stom_O3_ms, r_c, r_b)
%FSO3_MM Calculates the stomatal flux of ozone following the Mapping Manual (LRTAP 2017)
%
%   INPUTS:
%       O3_canopy_nmolm3  - Ozone concentration at canopy level [as given]
%       g_stom_O3_ms      - Stomatal conductance for ozone [m s^-1]
%       r_c               - Surface resistance of the entire leaf [s m^-1]
%       r_b               - Sublaminar resistance for a single leaf [s m^-1]
%
%   OUTPUT:
%       Fs_O3             - Instantaneous stomatal flux [same as input]
%
% Author: PR Guaita, 2025

Fs_O3 = O3 .* g_stom_O3_ms .* r_c ./ (r_b + r_c);

end
