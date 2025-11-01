function POD_Y = func_PODY(F_O3, Y, PODY_prev, timestep_s)
%FUNC_PODY Calculates the POD_Y metric (Phytotoxic Ozone Dose above a threshold Y)
%
%   INPUTS:
%       F_O3       - Instantaneous stomatal ozone flux [nmol m⁻² LAI s⁻¹]
%       Y          - Threshold flux [nmol m⁻² LAI s⁻¹]
%       PODY_prev  - Cumulative POD_Y up to previous timestep [mmol m⁻² LAI]
%       timestep_s - Duration of the time step [s]
%
%   OUTPUT:
%       POD_Y      - Updated cumulative POD_Y [mmol m⁻² LAI]
%
% Notes:
%   - Conversion from nmol to mmol accounts for unit consistency.
%   - Only fluxes above threshold Y contribute to the POD_Y.
%
% Author: PR Guaita, 2025

POD_Y = PODY_prev + max(0, (F_O3 - Y) * timestep_s / 1e6);

end
