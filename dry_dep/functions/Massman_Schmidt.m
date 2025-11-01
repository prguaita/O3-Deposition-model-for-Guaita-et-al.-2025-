function [ScO3, ScH2O, ScHeat] = Massman_Schmidt(T, P_z, config_const)
%MASSMAN_SCHMIDT Calculates Schmidt numbers for O3, H2O, and heat diffusion
%
%   INPUTS:
%     T          - Air temperature [°C]
%     P_z        - Atmospheric pressure [kPa]
%     config_const - Structure containing constants:
%                    D01_O3, D01_W, D01_Heat [m^2/s at reference conditions]
%                    alpha_diffusivity_O3, alpha_diffusivity_W, alpha_diffusivity_Heat [dimensionless exponents]
%
%   OUTPUTS:
%     ScO3       - Schmidt number for ozone (O3) diffusion [dimensionless]
%     ScH2O      - Schmidt number for water vapor (H2O) diffusion [dimensionless]
%     ScHeat     - Schmidt number for heat diffusion [dimensionless]
%
%   NOTES:
%     - Viscosity is estimated from temperature using Monteith (2014) data.
%     - The Schmidt number is calculated as: Sc = ν / (D * 100)
%       where ν is kinematic viscosity [cm^2/s] and D is molecular diffusivity [m^2/s].
%     - Temperature is converted to Kelvin for diffusivity calculations.
%
% Author: PR Guaita, 2025

D01_O3   = config_const.D01_O3;
D01_W    = config_const.D01_W;
D01_Heat = config_const.D01_Heat;
alpha_diffusivity_O3   = config_const.alpha_diffusivity_O3;
alpha_diffusivity_W    = config_const.alpha_diffusivity_W;
alpha_diffusivity_Heat = config_const.alpha_diffusivity_Heat;

T_0 = 273.16;  % Reference temperature [K]
P_0 = 101.325; % Reference pressure [kPa]

% Estimate kinematic viscosity (ν) in cm^2/s using Monteith (2014) data
v = zeros(size(T));
v_dummy = interp1(-5:5:45, [12.9 13.3 13.7 14.2 14.6 15.1 15.5 16 16.4 16.9 17.4], T);
v(T < -5 | T > 45) = 12.9 + T(T < -5 | T > 45) .* (17.4 - 12.9) / (45 - (-5));
v(T >= -5 & T <= 45) = v_dummy(T >= -5 & T <= 45);

% Calculate molecular diffusivities adjusted for temperature and pressure
D_O3   = D01_O3 .* (P_0 ./ P_z) .* (max(0, (T + 273.16) ./ T_0)).^(alpha_diffusivity_O3);
D_W    = D01_W .* (P_0 ./ P_z) .* (max(0, (T + 273.16) ./ T_0)).^(alpha_diffusivity_W);
D_heat = D01_Heat .* (P_0 ./ P_z) .* (max(0, (T + 273.16) ./ T_0)).^(alpha_diffusivity_Heat);

% Compute Schmidt numbers: Sc = kinematic viscosity / (diffusivity * 100)
ScO3   = v ./ (100 * D_O3);
ScH2O  = v ./ (100 * D_W);
ScHeat = v ./ (100 * D_heat);

end
