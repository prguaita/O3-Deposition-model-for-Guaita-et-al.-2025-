function [AWC, PAW] = bucket_model(AWC_prev, W_in_prev, E_prev, T_prev, AWHC, WaterSurplus)
%BUCKET_MODEL Computes soil water content recursively using a bucket model.
%
%   This function calculates the current Available Water Content (AWC) in the soil
%   using a simple bucket approach, accounting for inflows, outflows, and
%   water holding capacity. It also returns Plant Available Water (PAW) as a fraction.
%
% INPUTS:
%   AWC_prev     - Previous Available Water Content [mm]
%   W_in_prev    - Incoming water (e.g., rainfall not intercepted by canopy) [mm]
%   E_prev       - Evaporation loss [mm]
%   T_prev       - Transpiration loss [mm]
%   AWHC         - Available Water Holding Capacity of the soil [mm]
%   WaterSurplus - Additional water from root zone exploration or other sources [mm]
%
% OUTPUTS:
%   AWC - Updated Available Water Content [mm], limited to [0, AWHC]
%   PAW - Plant Available Water fraction (AWC divided by AWHC) [0-1]
%
% Author: PR Guaita, 2023

%% Update soil water content with water inputs and losses

AWC = AWC_prev + W_in_prev + WaterSurplus - E_prev - T_prev;

% Limit AWC to the physically possible range [0, AWHC]
AWC = max(0, min(AWC, AWHC));

%% Calculate Plant Available Water as a fraction

PAW = AWC ./ AWHC;

end
