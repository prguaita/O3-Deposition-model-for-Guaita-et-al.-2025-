function interp_map = interp_time(t, time_input, input_variable, option)
% interp_time Interpolates a 2D (or 1D) map at a specific time point.
%
%   interp_map = interp_time(t, time_input, input_variable, option)
%
%   Inputs:
%       t               - Target time at which to interpolate the map
%       time_input      - Vector of time coordinates corresponding to input_variable's third dimension
%       input_variable  - Data array (2D or 3D: lon x lat x time) of maps over time
%       option          - Interpolation method: 'linear', 'nearest', or 'mean'
%
%   Output:
%       interp_map      - Interpolated 2D map at time t
%
%   Description:
%       This function interpolates between two adjacent time slices in
%       input_variable to estimate the map at time t. If t is outside the
%       range of time_input, the closest available map is returned.
%
%   Author: PR Guaita - July 2023

    % Find the index of the first time in time_input greater than t
    ind_t_2 = find(t < time_input, 1, 'first');

    % If t is beyond the last available time, use the last time slice
    if isempty(ind_t_2)
        ind_t_2 = length(time_input);
    end

    % Handle the case when t is before or equal to the first time in time_input
    if ind_t_2 == 1
        % Extract the single map for interpolation (no interpolation needed)
        if ndims(input_variable) == 2
            interp_map = input_variable;
        else
            interp_map = input_variable(:, :, ind_t_2);
        end
        return
    end

    % Extract the two time points surrounding t
    t_1 = time_input(ind_t_2 - 1);
    t_2 = time_input(ind_t_2);

    % Extract the two maps corresponding to t_1 and t_2
    if ndims(input_variable) == 2
        Y_1 = input_variable;
        Y_2 = input_variable;
    else
        Y_1 = input_variable(:, :, ind_t_2 - 1);
        Y_2 = input_variable(:, :, ind_t_2);
    end

    % Perform interpolation based on the selected option
    switch option
        case 'linear'
            % Linear interpolation between Y_1 and Y_2
            interp_map = Y_1 + (t - t_1) * (Y_2 - Y_1) / (t_2 - t_1);
        case 'nearest'
            % Nearest neighbor interpolation: choose closer map
            if (t - t_1) < (t_2 - t)
                interp_map = Y_1;
            else
                interp_map = Y_2;
            end
        case 'mean'
            % Average of the two maps
            interp_map = (Y_1 + Y_2) / 2;
        otherwise
            error('Invalid option for interp_time. Choose ''linear'', ''nearest'', or ''mean''.');
    end

end
