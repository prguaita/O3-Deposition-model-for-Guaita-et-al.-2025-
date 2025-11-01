function [LAI_map, SAI_map, RD_map, Alb_map, h_c_map, ...
          d_map, z0M_map, z0H_map] = ...
    module_plant_phenology_v1(TT_data, whclim_map, config)
% MODULE_PLANT_PHENOLOGY_V1 Calculates plant parameters based on thermal time.
%
%   This function computes or loads plant parameters such as:
%   - LAI (Leaf Area Index)
%   - SAI (Surface Area Index)
%   - RD (Root Depth)
%   - Alb (Albedo)
%   - h_c (Canopy Height)
%   - d (Displacement height)
%   - z0M (Roughness length for momentum)
%   - z0H (Roughness length for heat)
%
%   The parameters are computed from thermal time data since sowing date,
%   with the option to use measured data or thermal-time-based models.
%
% Inputs:
%   TT_data     - 2D matrix of thermal time since sowing and up to 1st January [°C day]
%   whclim_map  - Climatic classification map to select parameter variants
%   config      - Configuration struct with fields:
%                   .plant.type        - 'measured' or 'TT'
%                   .plant.folder_plant- folder path for measured data
%                   .plant_TT          - thermal time parameters (winter)
%                   .plant_TT_spring   - thermal time parameters (spring)
%                   .plant_TT_mediterranean - thermal time params (mediterranean)
%                   .LAI, .SAI, .RD, .Alb - plant parameter structs or values
%                   .plant.h_c_max     - max canopy height
%                   .SAI.A_two         - normalization constant for SAI
%                   .plant.d_frac      - fraction to calculate displacement height
%                   .plant.z0M_frac    - fraction for momentum roughness length
%                   .plant.z0H_frac    - fraction for heat roughness length
%
% Outputs:
%   LAI_map     - Leaf Area Index map
%   SAI_map     - Stem Area Index map
%   RD_map      - Root Depth map
%   Alb_map     - Albedo map
%   h_c_map     - Canopy height map
%   d_map       - Displacement height map
%   z0M_map     - Roughness length for momentum
%   z0H_map     - Roughness length for heat

switch config.plant.type
    case 'measured'
        % Load plant parameters from stored files
        load(fullfile(main_folder, config.plant.folder_plant, 'LAI.mat'))
        load(fullfile(main_folder, config.plant.folder_plant, 'SAI.mat'))
        load(fullfile(main_folder, config.plant.folder_plant, 'RD.mat'))
        load(fullfile(main_folder, config.plant.folder_plant, 'Alb.mat'))
        load(fullfile(main_folder, config.plant.folder_plant, 'h_c.mat'))

    case 'TT'
        % Calculate plant parameters using thermal time models
        
        % Winter parameterization
        LAI_map = single(round(TT_based_var(TT_data, config.plant_TT, config.LAI, zeros(size(TT_data))), 2));
        RD_map  = single(round(RD_TT(TT_data, config.plant_TT, config.plant, zeros(size(TT_data))), 2));
        SAI_map = single(round(TT_based_var(TT_data, config.plant_TT, config.SAI, zeros(size(TT_data))), 2));
        Alb_map = single(round(config.plant.Alb * ones(size(TT_data)), 2));
        h_c_map = single(round(SAI_map * config.plant.h_c_max / config.SAI.A_two, 2));
        
        % Spring parameterization
        LAI_map_spring = single(round(TT_based_var(TT_data, config.plant_TT_spring, config.LAI, zeros(size(TT_data))), 2));
        RD_map_spring  = single(round(RD_TT(TT_data, config.plant_TT_spring, config.plant, zeros(size(TT_data))), 2));
        SAI_map_spring = single(round(TT_based_var(TT_data, config.plant_TT_spring, config.SAI, zeros(size(TT_data))), 2));
        h_c_map_spring = single(round(SAI_map_spring * config.plant.h_c_max / config.SAI.A_two, 2));
        
        % Mediterranean parameterization
        LAI_map_mediterranean = single(round(TT_based_var(TT_data, config.plant_TT_mediterranean, config.LAI, zeros(size(TT_data))), 2));
        RD_map_mediterranean  = single(round(RD_TT(TT_data, config.plant_TT_mediterranean, config.plant, zeros(size(TT_data))), 2));
        SAI_map_mediterranean = single(round(TT_based_var(TT_data, config.plant_TT_mediterranean, config.SAI, zeros(size(TT_data))), 2));
        h_c_map_mediterranean = single(round(SAI_map_mediterranean * config.plant.h_c_max / config.SAI.A_two, 2));
        
        % Substitute spring and mediterranean varieties in the corresponding climatic zones
        spring_idx = (whclim_map == 1.5);
        med_idx = (whclim_map == 0.5 | whclim_map == -0.5);
        
        LAI_map(spring_idx) = LAI_map_spring(spring_idx);
        RD_map(spring_idx) = RD_map_spring(spring_idx);
        SAI_map(spring_idx) = SAI_map_spring(spring_idx);
        h_c_map(spring_idx) = h_c_map_spring(spring_idx);
        
        LAI_map(med_idx) = LAI_map_mediterranean(med_idx);
        RD_map(med_idx) = RD_map_mediterranean(med_idx);
        SAI_map(med_idx) = SAI_map_mediterranean(med_idx);
        h_c_map(med_idx) = h_c_map_mediterranean(med_idx);
        
    otherwise
        error('Invalid config.plant.type. Must be either ''measured'' or ''TT''.');
end

% Calculate displacement and roughness length maps
d_map  = h_c_map * config.plant.d_frac;   % displacement height (d)
z0M_map = h_c_map * config.plant.z0M_frac; % roughness length for momentum (z0M)
z0H_map = z0M_map * config.plant.z0H_frac; % roughness length for heat (z0H)

disp('module Plant Phenology: DONE');

end
