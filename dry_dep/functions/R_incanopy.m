function R_inc = R_incanopy(SAI, u_star, h_c)
%R_INCANOPY Calculates the resistance inside the canopy [Jacob, 1994].
%
%   This function estimates the aerodynamic resistance within a vegetation canopy.
%
%   INPUTS:
%       SAI    [--]    Surface Area Index (dimensionless), matrix allowed
%       u_star [m/s]   Friction velocity, matrix allowed
%       h_c    [m]     Canopy height, matrix allowed
%
%   OUTPUT:
%       R_inc  [s/m]   Resistance inside the canopy, matrix output
%
%   Parameter:
%       b = 14 [m^-1], empirical coefficient from Jacob (1994)
%

b = 14; % Empirical coefficient [m^-1]

R_inc = b .* SAI .* h_c ./ u_star;

end
