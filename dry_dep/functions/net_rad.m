function R_n = net_rad(T_K, alpha, N, Alb, Q_sw)
%NET_RAD Computes net radiation at the surface.
%
%   Calculates the net radiation considering shortwave and longwave 
%   components and effects of temperature, albedo, and cloud cover.
%
%   INPUT:
%     T_K   - Air temperature [K]
%     alpha - Empirical parameter related to land cover / environment [unitless]
%     N     - Cloud cover fraction [unitless]
%     Alb   - Surface albedo [unitless]
%     Q_sw  - Incoming global shortwave radiation [W m^-2]
%
%   OUTPUT:
%     R_n   - Net radiation [W m^-2]
%
%   Notes:
%     - S depends exponentially on temperature and modifies the empirical 
%       correction factors (c3).
%     - sigma is the Stefan-Boltzmann constant.
%
% Author: PR Guaita, 2025

% Empirical parameters and physical constant
sigma = 5.67e-8;                % Stefan-Boltzmann constant [W m^-2 K^-4]
S = 1.5 * exp(-0.060208041 * T_K);
c1 = 5.31e-13; 
c2 = 60; 
c3 = 0.38 * ((1 - alpha) * S + 1) ./ (1 + S);

% Net radiation calculation
R_n = ((1 - Alb) .* Q_sw + c1 * T_K.^6 - sigma * T_K.^4 + c2 * N) ./ (1 + c3);

end
