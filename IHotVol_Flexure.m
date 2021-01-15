function [Xflx,Yflx,Zflx,HSPT_TRK]=IHotVol_Flexure(loadfile,rho_c,rho_w,rho_m,rho_i,T_elas,T_mant,kappa,HSPT_TRK,grdfile,ii)

%% Calculates a blended flexure profile using loading ages from HSPT_TRK

% generates blend file for blending flexure grids
system('rm Blendfile.txt');
fileID=fopen('Blendfile.txt','w');
formatspec='%s \n';
mkdir FLXcomptmp

% loading age
for jj=1:length(HSPT_TRK(:,1))

	Agediff=HSPT_TRK(jj,4)-HSPT_TRK(jj,3); % loading age = lith - seamount age

	if Agediff<0
   	Agediff=0.001; % litho min age is that of seamount (Te -> 0 if true)
	end

	Agediff=Agediff*1e6*365.25*24*3600; % unit conversion

	Te=erfinv(T_elas/T_mant)*2*sqrt(kappa*Agediff)/1e3; % elastic thickness

	HSPT_TRK(jj,6)=Te; % Te carried with track data

    % grd flexure calculation
	system(['grdflexure ' loadfile ' -D' num2str(rho_m) '/' ... 
    	num2str(rho_c) '/' num2str(rho_i) '/' num2str(rho_w) ...
    	' -E' num2str(Te) 'k -fg -N+a -GFLX_comp_Te_' ...
    	num2str(round(10*Te)/10) 'km.grd']);

	%append to blend file, with region
	Blendline = ['FLX_comp_Te_' num2str(round(10*Te)/10) 'km.grd -R' ...
    num2str(HSPT_TRK(jj,1)-2) '/' num2str(HSPT_TRK(jj,1)+2) '/' ...
    num2str(HSPT_TRK(jj,2)-2) '/' num2str(HSPT_TRK(jj,2)+2) ' ' num2str(1)];
  
	fprintf(fileID,formatspec,Blendline);

end

%close file
fclose(fileID);

%generate blended flexure grid
system(['grdblend Blendfile.txt -Gflexure.' num2str(ii) '.grd -R' grdfile ' -V']);

% stash pre-blend grids for later if needed
!mv FLX_comp* FLXcomptmp/

%DENAN flexure grid
system(['grdmath flexure.' num2str(ii) '.grd 0 DENAN = flexure.DENAN.' num2str(ii) '.grd']);

% read in flexure grid
[Xflx,Yflx,Zflx]=grdread2(['flexure.DENAN.' num2str(ii) '.grd']);

% --
