function [H, u_star, L, u_hc, ...
    R_aH_hc_zmO3, R_bO3, R_surf_O3, R_aH_dz0m_zmT, R_bW, R_bH, R_inc, ...
    r_s_act_O3, r_s_act_W, rc, rbO3, g_ext, ...
    O3_h_c_ugm3, O3_h_c_ppb, O3_2m_ppb, O3_3m_ppb, Fs_O3, Fs_O3_ppb, Ftot_O3_ppb, ...
    PODY, PODY_prev, POD0, POD0_prev, ...
    tleaf_C, VPD_leaf_kPa, e_sat_tleaf_kPa, ...
    tras_W, E_soil, E_wet, evap_W, evaptras_W, tras_W_prev, evap_W_prev, ...
    S_c, W_in, W_in_prev] = ...
    module_atmosphere_part_two(z_m_O3, tas_C, tas_K, e_sat_kPa, e_kPa, VPD_kPa, delta, ...
    u_z, rain, O3_zm_ugm3, H, ...
    light_flag, Rn_Wm2, Rn_soil, G_Wm2, N, ...
    P_z_hour, rho_air_kgm3, c_p_JkgK, Sc_W, Sc_H, Sc_O3, ...
    d, z_0m, LAI, SAI, h_c, PACC_map, ...
    g_s_act_O3_ms, g_s_act_W_ms, g_s_O3_sun_ms, PODY_prev, POD0_prev, ...
    AWC, S_c, W_in_prev, PAW_hour, ...
    config, par_g_stom)
% MODULE_ATMOSPHERE_PART_TWO
%
% Calculates ozone concentration at canopy height, aerodynamic and canopy resistances,
% stomatal conductance, friction velocity, and related fluxes.
%
% This module handles the atmospheric part of the surface energy and ozone exchange,
% including Monin–Obukhov stability, resistances, and deposition fluxes.
%
% Inputs:
%   z_m_O3         - Height of O3 measurement [m]
%   tas_C, tas_K   - Air temperature [°C, K]
%   e_sat_kPa      - Saturation vapor pressure [kPa]
%   e_kPa          - Actual vapor pressure [kPa]
%   VPD_kPa        - Vapor pressure deficit [kPa]
%   delta          - Slope of saturation vapor pressure curve [kPa K^-1]
%   u_z            - Wind speed at measurement height [m s^-1]
%   rain           - Rainfall [mm]
%   O3_zm_ugm3     - Ozone concentration at measurement height [µg m^-3]
%   H              - Sensible heat flux [W m^-2]
%   light_flag     - Boolean mask for grid cells with significant sunlight
%   Rn_Wm2         - Net radiation [W m^-2]
%   Rn_soil        - Net soil radiation [W m^-2]
%   G_Wm2          - Soil heat flux [W m^-2]
%   N              - Cloud cover fraction
%   P_z_hour       - Atmospheric pressure [hPa or kPa, depending on context]
%   rho_air_kgm3   - Air density [kg m^-3]
%   c_p_JkgK       - Specific heat capacity of air at constant pressure [J kg^-1 K^-1]
%   Sc_W, Sc_H, Sc_O3 - Schmidt numbers for water, heat, and ozone
%   d              - Displacement height [m]
%   z_0m           - Roughness length for momentum [m]
%   LAI            - Leaf Area Index
%   SAI            - Surface Area Index
%   h_c            - Canopy height [m]
%   PACC_map       - Plant accumulated canopy conductance map
%   g_s_act_O3_ms  - Actual stomatal conductance for ozone [m s^-1]
%   g_s_act_W_ms   - Actual stomatal conductance for water [m s^-1]
%   PODY_prev, POD0_prev - Previous timestep POD metrics
%   AWC            - Available water capacity
%   S_c            - Soil moisture stress coefficient
%   W_in_prev      - Previous water input
%   PAW_hour       - Plant available water fraction at current timestep
%   config         - Configuration struct with constants and options
%   par_g_stom     - Stomatal conductance parameters
%
% Outputs:
%   H, u_star, L            - Sensible heat flux, friction velocity, Monin–Obukhov length
%   u_hc                    - Wind speed at canopy height [m s^-1]
%   R_aH_hc_zmO3, R_bO3, R_surf_O3, R_aH_dz0m_zmT, R_bW, R_bH, R_inc - Various aerodynamic and boundary layer resistances [s m^-1]
%   r_s_act_O3, r_s_act_W  - Actual stomatal resistances [s m^-1]
%   rc, rbO3, g_ext        - Additional resistances and canopy conductance
%   O3_h_c_ugm3            - Ozone concentration at canopy height [µg m^-3]
%   O3_h_c_ppb, O3_2m_ppb, O3_3m_ppb - Ozone concentrations at canopy, 2 m and 3 m heights [ppb]
%   F_O3, Fs_O3_ppb, F_totDDIM_O3_ppb - Ozone fluxes and deposition metrics
%   PODY, POD0             - Updated phytotoxic ozone dose metrics
%   tleaf_C                - Estimated leaf temperature [°C]
%   VPD_leaf_kPa           - Vapor pressure deficit at leaf surface [kPa]
%   e_sat_tleaf_kPa        - Saturation vapor pressure at leaf temperature [kPa]
%   tras_W, E_soil, E_wet, evap_W, evaptras_W - Evapotranspiration and evaporation components
%   tras_W_prev, evap_W_prev - Previous timestep transpiration and evaporation
%   S_c, W_in, W_in_prev   - Updated soil moisture stress coefficient and water inputs
%
% Notes:
%   - Handles ozone scaling for different configurations.
%   - Includes Monin–Obukhov stability and resistances.
%   - Computes leaf-level processes and canopy fluxes.
%
% PR Guaita – 2025
    
%% --------------------------------------------------------------------------
% MODULE_ATMOSPHERE_PART_TWO - Part 1: Scale O3 concentration to canopy height
% Computes sensible heat flux, friction velocity, Monin–Obukhov length, 
% in-canopy resistance, and scales ozone concentration to canopy top.
%
% Author: PR Guaita, 2025
%--------------------------------------------------------------------------

[n_row, n_col] = size(PODY_prev);  % Grid dimensions

%% Step 1: Compute sensible heat flux (H)
% If 'neutral atmosphere' option is selected, assume H = 0
% Else, compute sensible heat flux if it is not provided as input
if config.dep.neutral_atm_option
    H = zeros(n_row, n_col);
else
    if ~config.input_model.H
        H = sens_heat(config.const, tas_K, light_flag, Rn_Wm2, G_Wm2);
    end
end

%% Step 2: Compute friction velocity (u_star) and Monin–Obukhov length (L)
% u_star depends on sensible heat flux, wind speed, canopy and atmospheric properties
u_star = fric_vel(H, u_z, d, z_0m, N, tas_K, rho_air_kgm3, c_p_JkgK, config.const);

% Check for numerical issues: complex u_star usually means negative wind speed somewhere
if ~isreal(u_star)
    disp('ERROR: u_star is complex. Likely due to negative wind speed in input data.');
    return
end

% Compute Monin–Obukhov length, L [m]
L = length_MO(rho_air_kgm3, tas_K, u_star, ...
    config.const.kar, config.const.g_0, H, c_p_JkgK);

%% Step 3: Compute in-canopy resistance
R_inc = R_incanopy(SAI, u_star, h_c);

%% Step 4: Scale ozone concentration from measurement height to canopy top
switch config.dep.ozone_scaling_option
    case 'measured'
        % Use measured O3 concentration directly at canopy height
        O3_h_c_ugm3 = O3_zm_ugm3;
        u_hc = u_z;

        % Other variables not computed in this option
        R_aH_hc_zmO3 = nan(n_row, n_col);
        R_bO3        = nan(n_row, n_col);
        R_surf_O3    = nan(n_row, n_col);
        r_s_act_O3   = nan(n_row, n_col);

    case 'same_tgt_ref'
        % Target (canopy top) and reference (measurement) heights differ
        % Compute wind speed at canopy top
        u_hc = scale_u_z(h_c, u_star, config.const.kar, d, z_0m, L);

        % Compute aerodynamic resistances between different heights
        R_aH_hc_zmO3   = R_aH(h_c, z_m_O3, d, config.const.kar, u_star, L);
        R_aH_2_zmO3    = R_aH(2, z_m_O3, d, config.const.kar, u_star, L);
        R_aH_3_zmO3    = R_aH(3, z_m_O3, d, config.const.kar, u_star, L);
        R_aH_hc_zmT    = R_aH(h_c, config.const.z_m_T, d, config.const.kar, u_star, L);
        R_aH_dz0m_zmO3 = R_aH(d+z_0m, z_m_O3, d, config.const.kar, u_star, L);

        % Compute boundary layer resistance for ozone
        R_bO3 = R_b(config.const, u_star, Sc_O3);

        % Stomatal resistance for ozone (inverse of conductance)
        r_s_act_O3 = 1 ./ g_s_act_O3_ms;

        % Surface resistance to ozone uptake (includes LAI, SAI, cuticular and soil resistances)
        R_surf_O3 = R_surf(LAI, SAI, r_s_act_O3, ...
            config.const.r_cut_O3, R_inc, config.const.R_soil);

        % Scale ozone concentration to canopy top and other heights
        O3_h_c_ugm3 = O3_z(O3_zm_ugm3, R_aH_hc_zmO3, ...
            R_aH_dz0m_zmO3 + R_bO3 + R_surf_O3);

        O3_2m_ugm3  = O3_z(O3_zm_ugm3, R_aH_2_zmO3, ...
            R_aH_dz0m_zmO3 + R_bO3 + R_surf_O3);

        O3_3m_ugm3  = O3_z(O3_zm_ugm3, R_aH_3_zmO3, ...
            R_aH_dz0m_zmO3 + R_bO3 + R_surf_O3);

    case 'different_tgt_ref'
        % Not implemented
        disp('ERROR: ozone_scaling_option = different_tgt_ref is not implemented.');
end

%--------------------------------------------------------------------------
% MODULE_ATMOSPHERE_PART_TWO - Part 2: Resistances, ozone fluxes, PODY, 
% evaporation and transpiration calculations
%
% Author: PR Guaita, 2025 (rewritten & documented)
%--------------------------------------------------------------------------

%% Step 5: Resistances (for heat, water and ozone)
% Aerodynamic resistance between roughness sublayer and measurement height (temperature)
R_aH_dz0m_zmT = R_aH(d+z_0m, config.const.z_m_T, d, ...
    config.const.kar, u_star, L);

% Boundary layer resistances (m s^-1) for water and heat
R_bW = R_b(config.const, u_star, Sc_W);
R_bH = R_b(config.const, u_star, Sc_H);

% Stomatal resistance for water (inverse of stomatal conductance)
r_s_act_W = 1 ./ g_s_act_W_ms;

% Surface resistance to ozone uptake by single leaf
g_ext = 1 / config.const.r_cut_O3;   % External cuticular conductance
rc = 1 ./ (g_s_O3_sun_ms + g_ext);  % Total surface resistance [s m^-1]

% Leaf boundary-layer resistances for ozone and heat transfer
rbO3 = 1.3 * 150 * sqrt(config.plant.cw_Leaf_dim ./ u_hc);  % O3
rbH  = 150 * sqrt(config.plant.cw_Leaf_dim ./ u_hc);        % Heat

%% Step 6: Convert O3 concentrations to ppb (from ug m^-3)
O3_h_c_ppb = ugm3_2_ppb(tas_K, config.const.M_O3, ...
    O3_h_c_ugm3, P_z_hour, config.const.Re);

O3_zm_ppb = ugm3_2_ppb(tas_K, config.const.M_O3, ...
    O3_zm_ugm3, P_z_hour, config.const.Re);

O3_2m_ppb = ugm3_2_ppb(tas_K, config.const.M_O3, ...
    O3_2m_ugm3, P_z_hour, config.const.Re);

O3_3m_ppb = ugm3_2_ppb(tas_K, config.const.M_O3, ...
    O3_3m_ugm3, P_z_hour, config.const.Re);

%% Step 7: Compute ozone fluxes
% Instantaneous stomatal ozone flux (ppb m s^-1) - following Mapping Manual
Fs_O3_ppb = FsO3_MM(O3_h_c_ppb, g_s_O3_sun_ms, rc, rbO3);

% Convert flux to nmol m^-2 s^-1
Fs_O3 = ppb_2_nmol_flux(Fs_O3_ppb, config.const.Re, P_z_hour * 1000, tas_K);

% Mask flux outside accumulation period
Fs_O3(~PACC_map) = 0;

%% Step 8: Compute PODY and POD0 (cumulative uptake)
PODY = func_PODY(Fs_O3, config.plant.Y, PODY_prev, config.timestep * 3600);  % mmol m^-2 LAI^-1
PODY_prev = PODY;

POD0 = func_PODY(Fs_O3, 0, POD0_prev, config.timestep * 3600);
POD0_prev = POD0;

% Total ozone flux at canopy scale (DDIM framework)
Ftot_O3_ppb = FtotO3_DDIM(O3_zm_ppb, R_aH_dz0m_zmO3, R_bO3, R_surf_O3);

%% Step 9: Stomatal transpiration
tras_W = tras_stom_mms(tas_C, rho_air_kgm3, c_p_JkgK, ...
    e_sat_kPa - e_kPa, R_aH_dz0m_zmT, R_bH, delta, Rn_Wm2, G_Wm2, ...
    config.plant.n, config.const.gamma, R_bW, r_s_act_W, LAI);

% Fix numeric issue
tras_W(isnan(tras_W)) = 0;

% Scale by timestep and soil water content
if isnan(AWC)
    tras_W = config.timestep_s * tras_W;
else
    tras_W = min(max(AWC - W_in_prev, 0), config.timestep * 3600 * tras_W);
end

%% Step 10: Soil evaporation
switch par_g_stom.soil_option
    case {'PAW', 'SWC'}
        if PAW_hour < 0.5
            E_soil = zeros(n_row, n_col);
        else
            E_soil = evap_soil_mms(tas_C, rho_air_kgm3, c_p_JkgK, ...
                VPD_kPa, R_aH_dz0m_zmT, R_bH, R_inc, ...
                delta, Rn_soil, G_Wm2, 2, config.const.gamma, R_bW, config.const.R_soil);
            if isnan(AWC)
                E_soil = config.timestep_s * E_soil;
            else
                E_soil = min(max(AWC - W_in_prev, 0), config.timestep * 3600 * E_soil);
            end
        end
    case 'FC'
        E_soil = zeros(n_row, n_col);
end

%% Step 11: Evaporation from wet surfaces (e.g. leaf wetness)
E_wet = min(S_c, config.timestep * 3600 * ...
    evap_sur_wet_mms(tas_C, rho_air_kgm3, c_p_JkgK, ...
        VPD_kPa, R_aH_dz0m_zmT, R_bH, delta, Rn_Wm2, G_Wm2, ...
        2, config.const.gamma, R_bW));

% Update water stored on leaves (S_c) and water infiltrated into soil (W_in)
[S_c, W_in] = rain_on_leaves(LAI, rain, S_c, E_wet);

%% Step 12: Aggregate total evaporation and evapotranspiration
evap_W     = E_soil + E_wet;    % Total evaporation
evaptras_W = evap_W + tras_W;   % Total evapotranspiration

%% Step 13: Store values depending on accumulation option
if config.dep.ET_ACC_option
    tras_W_prev = zeros(n_row, n_col);
    evap_W_prev = zeros(n_row, n_col);
    tras_W_prev(PACC_map) = tras_W(PACC_map);
    evap_W_prev(PACC_map) = evap_W(PACC_map);
else
    tras_W_prev = tras_W;
    evap_W_prev = evap_W;
end

%% Step 14: Sanity check to avoid complex numbers
if ~isreal(POD0) || ~isreal(PODY) || ~isreal(Fs_O3) ...
   || ~isreal(g_s_act_O3_ms) || ~isreal(PAW_hour) || ~isreal(tras_W)
    disp(['Complex values detected at hour ' int2str(hoy)]);
    return
end

% Update previous water input for next timestep
W_in_prev = W_in;

end