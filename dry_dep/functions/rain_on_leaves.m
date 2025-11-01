function [S_c, W_in] = rain_on_leaves(LAI, rain, S_c_prev, E_wet)
%RAIN_ON_LEAVES Computes canopy interception and throughfall.
%
%   This function calculates how much rainwater is stored on the leaf surface 
%   and how much directly reaches the soil after exceeding storage capacity.
%
%   INPUT:
%     LAI       - Leaf Area Index [m^2 leaf area per m^2 ground area]
%     rain      - Rainfall during the current timestep [mm]
%     S_c_prev  - Previous canopy water storage [mm]
%     E_wet     - Evaporation from wet surfaces during the timestep [mm]
%
%   OUTPUT:
%     S_c       - Updated canopy water storage [mm]
%     W_in      - Water input to the soil (throughfall and drainage) [mm]
%
%   Notes:
%     - Canopy maximum storage is assumed as 0.1 mm per unit of LAI.
%     - Negative storage or water input is forced to zero.
%
% Author: PR Guaita, 2025

% Maximum canopy storage capacity [mm]
S_c_max = 0.1 * LAI;

% Updated canopy storage: cannot exceed capacity, cannot be negative
S_c = min(S_c_max, max(rain + S_c_prev - E_wet, 0));

% Water input to soil: rain exceeding storage capacity
W_in = max(rain + S_c_prev - S_c_max, 0);

end
