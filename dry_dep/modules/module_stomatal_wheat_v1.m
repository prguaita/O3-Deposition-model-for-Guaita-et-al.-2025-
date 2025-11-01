function [g_s_frac_sun, g_s_frac_shad, sumVPD, jarvis_func, Uddling, ...
          g_s_W_sun_ms, g_s_O3_sun_ms, g_s_W_shad_ms, g_s_O3_shad_ms, ...
          g_s_act_W_ms, g_s_act_O3_ms, g_s_frac_sun_prev, g_s_frac_shad_prev] = ...
          module_stomatal_wheat_v1(tas, VPD, PAW, PPFD_sun, PPFD_shad, CO2_ppm, ...
          POD0_prev, P_z_hour, sumVPD, ...
          g_s_frac_sun_prev, g_s_frac_shad_prev, g_max_W_winter, g_max_O3_winter, ...
          g_max_W_medite, g_max_O3_medite, LAI_sun, LAI_shad, LAI, TT_data, ...
          par_g_stom_winter, par_g_stom_spring, par_g_stom_medite, whclim_map, config)
% MODULE_STOMATAL_WHEAT_V1 Calculates stomatal conductance fractions and scaled conductances for wheat.
%
%   Outputs stomatal conductance fractions for sunlit and shaded leaves,
%   scaled stomatal conductances for water vapor and ozone, and updates
%   internal state variables for use in subsequent time steps.
%
% Inputs:
%   tas             - Air temperature (°C)
%   VPD             - Vapor Pressure Deficit (kPa)
%   PAW             - Plant Available Water
%   PPFD_sun        - Photosynthetic Photon Flux Density for sunlit leaves
%   PPFD_shad       - PPFD for shaded leaves
%   CO2_ppm         - CO2 concentration (ppm)
%   POD0_prev       - Previous ozone phytotoxic dose
%   P_z_hour        - Pressure at height (Pa)
%   sumVPD          - Cumulative VPD from prior timesteps
%   g_s_frac_sun_prev - Previous stomatal conductance fraction for sunlit leaves
%   g_s_frac_shad_prev - Previous stomatal conductance fraction for shaded leaves
%   g_max_W_winter  - Max stomatal conductance for water vapor (winter params)
%   g_max_O3_winter - Max stomatal conductance for ozone (winter params)
%   g_max_W_medite  - Max stomatal conductance for water vapor (Mediterranean params)
%   g_max_O3_medite - Max stomatal conductance for ozone (Mediterranean params)
%   LAI_sun         - Leaf Area Index sunlit leaves
%   LAI_shad        - Leaf Area Index shaded leaves
%   LAI             - Total Leaf Area Index
%   TT_data         - Thermal time map (2D)
%   par_g_stom_winter, par_g_stom_spring, par_g_stom_medite - Parameter structs for stomatal conductance
%   whclim_map      - Climate classification map (used to flag Mediterranean nodes)
%   config          - Configuration struct
%
% Outputs:
%   g_s_frac_sun, g_s_frac_shad - Fraction of max stomatal conductance (sunlit/shaded)
%   sumVPD                      - Updated cumulative VPD
%   jarvis_func, Uddling        - Structs with intermediate stomatal function outputs
%   g_s_W_sun_ms, g_s_O3_sun_ms - Stomatal conductance for water vapor and ozone [m/s] (sunlit)
%   g_s_W_shad_ms, g_s_O3_shad_ms - Same as above but for shaded leaves
%   g_s_act_W_ms, g_s_act_O3_ms - LAI-weighted active stomatal conductance for water and ozone
%   g_s_frac_sun_prev, g_s_frac_shad_prev - Updated fraction values for next timestep

    % Calculate stomatal conductance fractions and related variables
    [g_s_frac_sun, g_s_frac_shad, sumVPD, jarvis_func, Uddling] = ...
        conductance_wheat_v1(TT_data, par_g_stom_winter, par_g_stom_spring, par_g_stom_medite, ...
        whclim_map, config, POD0_prev, PPFD_sun, PPFD_shad, tas, VPD, PAW, CO2_ppm, ...
        sumVPD, g_s_frac_sun_prev, g_s_frac_shad_prev);

    % Update previous fractions for next timestep
    g_s_frac_sun_prev = g_s_frac_sun;
    g_s_frac_shad_prev = g_s_frac_shad;

    % Calculate stomatal conductance for H2O and O3 (winter parameters)
    g_s_W_sun  = g_max_W_winter  .* g_s_frac_sun;  % mmol H2O m^-2 PLA s^-1
    g_s_O3_sun = g_max_O3_winter .* g_s_frac_sun;  % mmol O3 m^-2 PLA s^-1
    g_s_W_shad = g_max_W_winter  .* g_s_frac_shad; % mmol H2O m^-2 PLA s^-1
    g_s_O3_shad= g_max_O3_winter .* g_s_frac_shad; % mmol O3 m^-2 PLA s^-1

    % Calculate stomatal conductance for Mediterranean parametrization
    g_s_W_sun_tmp  = g_max_W_medite  .* g_s_frac_sun;
    g_s_O3_sun_tmp = g_max_O3_medite .* g_s_frac_sun;
    g_s_W_shad_tmp = g_max_W_medite  .* g_s_frac_shad;
    g_s_O3_shad_tmp= g_max_O3_medite .* g_s_frac_shad;

    % Flag Mediterranean climate nodes
    flag_medite_wheat = (whclim_map == -0.5) | (whclim_map == 0.5);

    % Replace winter gs with Mediterranean gs for flagged nodes
    g_s_W_sun(flag_medite_wheat)   = g_s_W_sun_tmp(flag_medite_wheat);
    g_s_O3_sun(flag_medite_wheat)  = g_s_O3_sun_tmp(flag_medite_wheat);
    g_s_W_shad(flag_medite_wheat)  = g_s_W_shad_tmp(flag_medite_wheat);
    g_s_O3_shad(flag_medite_wheat) = g_s_O3_shad_tmp(flag_medite_wheat);

    % Convert stomatal conductance from mmol m^-2 s^-1 to m/s
    g_s_W_sun_ms  = mmol2ms(g_s_W_sun, P_z_hour, tas);
    g_s_O3_sun_ms = mmol2ms(g_s_O3_sun, P_z_hour, tas);
    g_s_W_shad_ms = mmol2ms(g_s_W_shad, P_z_hour, tas);
    g_s_O3_shad_ms= mmol2ms(g_s_O3_shad, P_z_hour, tas);

    % Calculate LAI-weighted stomatal conductance (water and ozone)
    g_s_act_W_ms  = (LAI_sun .* g_s_W_sun_ms + LAI_shad .* g_s_W_shad_ms) ./ LAI;
    g_s_act_O3_ms = (LAI_sun .* g_s_O3_sun_ms + LAI_shad .* g_s_O3_shad_ms) ./ LAI;

end
