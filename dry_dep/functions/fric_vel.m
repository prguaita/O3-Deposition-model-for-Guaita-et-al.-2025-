function u_star = fric_vel(H, u_z, d, z_0, N, T_K, rho, c_p, config_const)
%FRIC_VEL Calculates the friction velocity (u*) following Hanna and Chang (1989).
%
%   Friction velocity (u*) is a fundamental parameter in boundary layer meteorology,
%   representing the shear velocity at the surface, related to momentum transfer.
%
%   INPUTS:
%       H           [W/m^2]   Sensible heat flux (positive = unstable, negative = stable)
%       u_z         [m/s]     Wind speed measured at height z
%       d           [m]       Zero-plane displacement height (canopy or surface displacement)
%       z_0         [m]       Roughness length of the surface
%       N           [--]      Cloud cover fraction (dimensionless, 0 to 1)
%       T_K         [K]       Air temperature in Kelvin
%       rho         [kg/m^3]  Air density
%       c_p         [J/(kg·K)] Specific heat capacity of air at constant pressure
%       config_const struct     Struct containing constants:
%                   - z_m_wnd: measurement height of wind speed (m)
%                   - kar: von Karman constant (~0.4)
%                   - g_0: acceleration due to gravity (m/s^2)
%
%   OUTPUT:
%       u_star      [m/s]     Friction velocity
%
% PR Guaita – 2023

% Extract constants from config struct
z_m_wnd = config_const.z_m_wnd; % Wind measurement height [m]
k       = config_const.kar;     % von Karman constant (~0.4)
g       = config_const.g_0;     % Gravity acceleration [m/s^2]

% Identify atmospheric stability based on sensible heat flux H
stable_mask   = H < -1;            % Stable atmosphere
neutral_mask  = H >= -1 & H < 1;   % Neutral atmosphere
unstable_mask = H >= 1;             % Unstable atmosphere

% Initialize output friction velocity vector
u_star = zeros(size(u_z));

% Calculate the logarithmic wind profile term
log_term = log((z_m_wnd - d) ./ z_0);

%% Stable atmospheric conditions
% Calculate the stability parameter theta_star with upper limit influenced by cloud cover
theta_star = min(0.09 * (1 - 0.5 * N.^2), ...
                 (k * T_K .* (u_z.^2)) ./ (18.8 * g * z_m_wnd .* log_term));

% Intermediate terms for stable friction velocity calculation
f1 = (0.5 * k * u_z) ./ log_term;
f2 = 1 - 4 * (4.7 * g * z_m_wnd .* theta_star .* log_term) ./ (k * T_K .* (u_z.^2));

% Address numerical precision: clamp values near zero to zero
f2(abs(f2) < 1e-12) = 0;

% Calculate friction velocity for stable conditions
u_star(stable_mask) = f1(stable_mask) .* (1 + sqrt(f2(stable_mask)));

%% Neutral atmospheric conditions
% Classic log wind profile formula for neutral conditions
u_star(neutral_mask) = (k * u_z(neutral_mask)) ./ log_term(neutral_mask);

%% Unstable atmospheric conditions
% Correction factors for unstable atmospheric conditions
d_f = z_0 ./ (z_m_wnd - d);
d1 = zeros(size(u_z));

% Assign empirical d1 values based on roughness fraction
d1(d_f <= 0.01) = 0.128 + 0.005 ./ log_term(d_f <= 0.01);
d1(d_f > 0.01) = 0.107;

% Calculate d2 parameter empirically based on roughness fraction
d2 = 1.95 + 32.6 * d_f.^0.45;

% Calculate d3 term incorporating heat flux and meteorological parameters
d3 = (H ./ (rho .* c_p)) .* (k * g * (z_m_wnd - d) ./ T_K) .* (log_term ./ (k * u_z)).^3;

% Compute friction velocity for unstable conditions
u_star(unstable_mask) = (k * u_z(unstable_mask) ./ log_term(unstable_mask)) .* ...
    (1 + d1(unstable_mask) .* log(1 + d2(unstable_mask) .* d3(unstable_mask)));

end
