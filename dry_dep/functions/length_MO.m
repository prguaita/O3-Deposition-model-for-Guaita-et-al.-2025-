function L = length_MO(rho, T_K, u_star, kar, g, H, c_p)
%LENGTH_MO Calculates the Monin-Obukhov length (L).
%
%   The Monin-Obukhov length is a fundamental parameter in micrometeorology
%   that characterizes the stability of the atmospheric surface layer.
%
%   INPUTS:
%       rho     [kg/m^3]    Air density
%       T_K     [K]         Air temperature in Kelvin
%       u_star  [m/s]       Friction velocity
%       kar     [--]        von Karman constant (~0.4)
%       g       [m/s^2]     Acceleration due to gravity
%       H       [W/m^2]     Sensible heat flux
%       c_p     [J/(kg·K)]  Specific heat capacity of air at constant pressure
%
%   OUTPUT:
%       L       [m]         Monin-Obukhov length
%

L = - (u_star).^3 ./ ( (kar * g * H) ./ (rho .* c_p .* T_K) );

end
