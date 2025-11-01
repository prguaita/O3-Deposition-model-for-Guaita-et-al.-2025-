function [PAR_tot_act, PAR_dir, PAR_diff, ...
          PAR_sun, PAR_shad, light_flag, ...
          PPFD_sun, PPFD_shad] = ...
    module_light_v1(B, Q_sw, P_z_hour, LAI, ...
                   nRow, nCol, config)
% MODULE_LIGHT_V1 Computes photosynthetically active radiation (PAR) components.
%
%   This function calculates the total, direct, and diffuse PAR, as well as
%   the sunlit and shaded leaf PAR components, and corresponding PPFD values.
%   It also generates a light flag indicating sufficient light levels.
%
% Inputs:
%   B           - Solar elevation (°)
%   Q_sw        - Incoming shortwave radiation [W m^-2]
%   P_z_hour    - Solar zenith angle or related parameter for PAR partitioning
%   LAI         - Leaf Area Index (dimensionless)
%   nRow, nCol  - Dimensions of the spatial grid (rows and columns)
%   config      - Configuration struct containing constants, including:
%                   .const.Qsw2PAR  - Conversion factor from shortwave to PAR
%                   .const.Wm2_PPFD - Conversion factor from W m^-2 to PPFD (µmol m^-2 s^-1)
%
% Outputs:
%   PAR_tot_act - Total photosynthetically active radiation [W m^-2]
%   PAR_dir     - Direct component of PAR [W m^-2]
%   PAR_diff    - Diffuse component of PAR [W m^-2]
%   PAR_sun     - PAR absorbed by sunlit leaves [W m^-2]
%   PAR_shad    - PAR absorbed by shaded leaves [W m^-2]
%   light_flag  - Logical array indicating light presence (true if > 50 W m^-2 and B>0)
%   PPFD_sun    - PPFD absorbed by sunlit leaves [µmol m^-2 s^-1]
%   PPFD_shad   - PPFD absorbed by shaded leaves [µmol m^-2 s^-1]

% Initialize light flag: true where Q_sw > 50 W m^-2 and solar elevation > 0
light_flag = false(nRow, nCol);
light_flag(Q_sw > 50 & B > 0) = true;

% Convert shortwave radiation to PAR
PAR_tot_act = config.const.Qsw2PAR * Q_sw;

% Calculate direct and diffuse PAR fractions
[PAR_dir, PAR_diff] = PAR_fractions(PAR_tot_act, P_z_hour, B);

% Calculate PAR absorbed by sunlit and shaded leaves
[PAR_sun, PAR_shad] = PAR_leaf_fractions(PAR_dir, PAR_diff, LAI, B, 60);

% Set PAR to zero where there is no sufficient light
PAR_sun(~light_flag) = 0;
PAR_shad(~light_flag) = 0;

% Convert PAR (W m^-2) to PPFD (µmol m^-2 s^-1)
PPFD_sun  = PAR_sun * config.const.Wm2_PPFD;
PPFD_shad = PAR_shad * config.const.Wm2_PPFD;

end
