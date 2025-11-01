function RD_map = RD_TT(TT,config_TT,config_plant,RD_map)
%RD_TT calculates the Root depth based on accumulated thermal time 
%   INPUT:
%   TT          - thermal Time from sowing
%   config      - struct variable that contains information about LAI
%   RD_Map      - preset Root depth map

TT_A_end   = config_TT.TT_A_end;

ndoy = size(RD_map,3);
for i_doy=1:ndoy
    HUI=TT(:,:,i_doy)./TT_A_end;
    RD_map(:,:,i_doy)=min(1,2.5*config_plant.RD_max*HUI);
end

end

