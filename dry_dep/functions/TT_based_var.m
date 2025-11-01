function var_map = TT_based_var(TT,config_TT,config_var,var_map)
%TT_based_var calculates a variable based on thermal time and specified
%parameters
%   INPUT:
%   TT          - thermal Time from sowing
%   config      - struct variable that contains information about LAI
%   var_map     - preset variable map


%% TT patameters

TT_one     = config_TT.TT_one; 
TT_two     = config_TT.TT_two; 
TT_three   = config_TT.TT_three;
TT_A_start = config_TT.TT_A_start;
TT_A_two   = config_TT.TT_A_two;
TT_A_three = config_TT.TT_A_three;
TT_A_end   = config_TT.TT_A_end;

%corresponding values
y1     = config_var.one; 
y2     = config_var.two; 
y3     = config_var.three; 
y_Astart = config_var.A_start;
y_A2   = config_var.A_two;
y_A3 = config_var.A_three;
y_Aend   = config_var.A_end;

TT_array    = [ 0 , TT_one , TT_two , TT_three , TT_A_start , TT_A_two , TT_A_three , TT_A_end  , TT_A_end+1 , 100000];
var_array   = [ 0 , y1, y2, y3  , y_Astart, y_A2, y_A3, y_Aend , 0          , 0] ;


%% Calcolo variabile

ndoy = size(var_map,3);
for doy=1:ndoy
	% interpolate
	var_map(:,:,doy)=interp1(TT_array,var_array,TT(:,:,doy));
end

end
