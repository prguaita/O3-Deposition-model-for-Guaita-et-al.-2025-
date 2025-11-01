function conc_ugm3 = ppb_2_ugm3(T, peso_mol, conc_ppb, P_z)
%PPB_2_UGM3 Converts gas concentration from ppb to micrograms per cubic meter [µg/m³]
%
%   INPUTS:
%     T         - Temperature in degrees Celsius [°C]
%     peso_mol  - Molecular weight of the gas [g/mol]
%     conc_ppb  - Gas concentration in parts per billion [ppb]
%     P_z       - Atmospheric pressure [kPa]
%
%   OUTPUT:
%     conc_ugm3 - Gas concentration in micrograms per cubic meter [µg/m³]
%
%   Note: Uses the ideal gas law for conversion.

conc_ugm3 = conc_ppb ./ ( (8.314 * (T + 273.15)) * 1000 ./ (peso_mol .* (P_z * 1000)) );

end
