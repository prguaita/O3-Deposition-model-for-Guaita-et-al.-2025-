function SAI_map = SAI_TT(TT,config,SAI_map)
%SAI_TT calculates the SAI based on accumulated thermal time
%   INPUT:
%   TT          - thermal Time from sowing
%   config      - struct variable that contains information about SAI
%   SAI_map     - preset SAI map


%% Definizione parametri

% TT semina-raccolto per le fasi del SAI
TT_one     = config.TT_one; 
TT_two     = config.TT_two; 
TT_three   = config.TT_three;
TT_A_start = config.TT_A_start;
TT_A_two   = config.TT_A_two;
TT_A_three = config.TT_A_three;
TT_A_end   = config.TT_A_end;

%valori del SAI corrispondenti
SAI_uno     = config.SAI_uno; 
SAI_due     = config.SAI_due; 
SAI_tre     = config.SAI_tre; 
SAI_A_start = config.SAI_A_start;
SAI_A_two   = config.SAI_A_two;
SAI_A_three = config.SAI_A_three;
SAI_A_end   = config.SAI_A_end;

TT_array    = [ 0 , TT_one , TT_two , TT_three , TT_A_start , TT_A_two , TT_A_three , TT_A_end  , TT_A_end+1, 100000];
SAI_array   = [ 0 , SAI_uno, SAI_due, SAI_tre  , SAI_A_start, SAI_A_two, SAI_A_three, SAI_A_end , 0         , 0] ;


%% Calcolo SAI
for doy=1:ndoy
	% calcola il SAI con un interpolazione di primo grado
	SAI_map(:,:,doy)=interp1(TT_array,SAI_array,TT(:,:,doy));
end

end

