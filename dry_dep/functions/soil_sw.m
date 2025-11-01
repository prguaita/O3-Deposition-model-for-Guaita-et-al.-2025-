function Q_sw_soil = soil_sw(Q_sw, Ka, LAI)
%SOIL_SW Computes the shortwave radiation reaching the soil surface.
%
%   Calculates the fraction of incoming global shortwave radiation (Q_sw) 
%   that penetrates the canopy and reaches the soil, based on the Beer-Lambert 
%   law and canopy attenuation.
%
%   INPUT:
%     Q_sw - Incoming global shortwave radiation above canopy [W m^-2]
%     Ka   - Canopy extinction coefficient for shortwave radiation [unitless]
%     LAI  - Leaf Area Index [m^2 leaf m^-2 ground]
%
%   OUTPUT:
%     Q_sw_soil - Shortwave radiation reaching the soil surface [W m^-2]
%
%   Formula:
%     Q_sw_soil = Q_sw * exp(-Ka * LAI)
%
% Author: PR Guaita, 2025

Q_sw_soil = Q_sw .* exp(-Ka * LAI);

end
