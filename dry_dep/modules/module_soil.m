function [W_cont, SWC_hour, PAW_hour, AWC_prev, AWC] = ...
    module_soil(AWC_prev, W_in_prev, evap_W_prev, tras_W_prev, root_surplus, ...
                soil_water_index, MWHC, config, par_g_stom, WP, FC)
% MODULE_SOIL Calculates soil water content and available water capacity.
%
%   This module calculates Plant Available Water (PAW), Available Water Capacity (AWC),
%   and Soil Water Content (SWC) based on the soil option selected in the configuration.
%
% Inputs:
%   AWC_prev        - Available Water Capacity at previous timestep
%   W_in_prev       - Previous water input to soil (e.g., rainfall, irrigation)
%   evap_W_prev     - Previous evaporation water loss
%   tras_W_prev     - Previous transpiration water loss
%   root_surplus    - Root zone surplus water
%   soil_water_index- Soil water index (measured or modeled)
%   MWHC            - Maximum Water Holding Capacity of the soil
%   config          - Configuration struct with .dep.soil_option and other options
%   par_g_stom      - Parameters related to stomatal conductance, including soil option
%   WP              - Wilting Point (soil water content)
%   FC              - Field Capacity (soil water content)
%
% Outputs:
%   W_cont          - Soil water content used for conductance calculation (SWC or PAW)
%   SWC_hour        - Soil Water Content (dimensionless fraction)
%   PAW_hour        - Plant Available Water (dimensionless fraction)
%   AWC_prev        - Updated Available Water Capacity (for next timestep)
%   AWC             - Current Available Water Capacity

switch config.dep.soil_option
    case 'full'
        % Full soil data: set AWC to max water holding capacity
        AWC = MWHC;
        PAW_hour = ones(size(AWC)); % Assume full PAW
        SWC_hour = ones(size(AWC)); % Assume full SWC
        W_cont = PAW_hour;          % Default water content for conductance

        % Choose water content for conductance based on stomatal soil option
        switch par_g_stom.soil_option
            case 'SWC'
                W_cont = SWC_hour;
            case 'PAW'
                W_cont = PAW_hour;
            otherwise
                % Optional: handle unexpected options here
        end

    case 'measured'
        % Use measured soil water index data
        if config.data.SWC
            SWC_hour = soil_water_index;
            AWC = nan;
            PAW_hour = nan;
            W_cont = SWC_hour;
        elseif config.data.PAW
            SWC_hour = nan;
            AWC = nan;
            PAW_hour = soil_water_index;
            W_cont = PAW_hour;
        else
            error('Measured soil data requires either SWC or PAW to be enabled in config.data');
        end

    case 'model'
        % Model soil water content using water balance model
        [AWC, PAW_hour] = bucket_model(AWC_prev, W_in_prev, evap_W_prev, tras_W_prev, MWHC, root_surplus);
        SWC_hour = AWC2SWC(MWHC, AWC, WP, FC);

        switch par_g_stom.soil_option
            case 'PAW'
                W_cont = PAW_hour;
            otherwise
                error('ERROR: only PAW can be provided for modeled soil data');
        end

    otherwise
        error('Unknown soil_option in config.dep.soil_option');
end

% Update AWC_prev for next timestep
AWC_prev = AWC;

end
