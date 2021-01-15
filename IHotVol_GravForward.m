function [XgMod YgMod ZgMod]=IHotVol_GravForward(grdfile,denangrdfile,edificegrdfile,sedcutgrdfile,subaerialgrdfile,mohoflexgrdfile,rho_c,rho_w,rho_m,rho_i,kappa,INP,ORS_L,subaq)

%% Forward gravity model

% gravity contribution from sed/water interface
system(['gravfft ' denangrdfile ' -D' num2str(rho_i-rho_w) ' -E5 -N+a -fg -Ggravmodel/sed.grav.grd']);
        
% gravity contribution from rock/sed & rock/air interface
if subaq==0
    
    % variable density grid for subaerial/subaqueous portions
    system(['grdmath ' subaerialgrdfile ' ' num2str(2400) ' MUL = subaerial_rho.grd']);
    system(['grdmath ' subaerialgrdfile ' 1 SUB -1 DIV ' num2str(rho_c-rho_i) ' MUL subaerial_rho.grd ADD = rho_var.grd']);

	% gravity from interface    
    system(['gravfft ' denangrdfile ' -Drho_var.grd -E4 -N+a -fg -Ggravmodel/rock.grav.grd']);
    system(['grdmath ' denangrdfile ' ' subaerialgrdfile ' MUL = subaerialtopo.grd']);
    system(['grdmath subaerialtopo.grd 0 LE gravmodel/sed.grav.grd MUL = gravmodel/sed.grav.grd']);
          
end

if subaq==1
    system(['grdmath ' denangrdfile ' ' sedcutgrdfile ' SUB = gravmodel/rockoffset.grd']);
    system(['gravfft gravmodel/rockoffset.grd -D' num2str(rho_c-rho_i) ' -E4 -N+a -fg -Ggravmodel/rock.grav.grd']);
end

% moho gravity (bottom of 7km crust, 5km of water, plus flexure)
system(['grdmath ' mohoflexgrdfile ' 12000 SUB = ' mohoflexgrdfile ]);
system(['gravfft ' mohoflexgrdfile ' -D' ...
    num2str(rho_m-rho_c) ' -E4 -fg -N+a -Ggravmodel/flex.moho.grav.grd']);


% get age and regional bathymetry
system(['grdsample cutagetmp/age.' grdfile ' -R' grdfile '_regional.grd -Ggravmodel/age.trim.' grdfile]);

% gravity contribution from 50C slices up to 600 degree isotherm
	for ii=50:50:550
    	%z=erfinv(ii/600)*2*sqrt(kappa)*sqrt(t);
    	COEF=erfinv(ii/600)*2*sqrt(kappa);
    	delRho=(50)*3e-5*3300;
    	ageconv=1e6*365.25*24*3600;

    	system(['grdmath ' grdfile '_regional.grd gravmodel/age.trim.' grdfile ' '...
       		num2str(ageconv) ' MUL SQRT ' ...
       		num2str(COEF) ' MUL SUB = gravmodel/depth.' num2str(ii) 'C.grd']);
    	system(['grdmath gravmodel/depth.' num2str(ii) 'C.grd 0 DENAN = gravmodel/depth.' num2str(ii) 'C.grd']);

    	if ii==50
			
		    !rm thermal.grav.grd
		    
		   	system(['gravfft gravmodel/depth.' num2str(ii) 'C.grd -D' num2str(-delRho) ...
		   		' -E4 -fg -N+a -Ggravmodel/thermal.grav.grd']);
		else
		   	system(['gravfft gravmodel/depth.' num2str(ii) 'C.grd -D' num2str(-delRho) ...
		   		' -E4 -fg -N+a -Ggravmodel/Grav.tmp.grd']);
		   	system('grdmath gravmodel/Grav.tmp.grd gravmodel/thermal.grav.grd ADD = gravmodel/thermal.grav.grd');
		   	
    	end
    
	end

% resize grids
system('grdsample -Rgravmodel/thermal.grav.grd gravmodel/flex.moho.grav.grd -Ggravmodel/flex.moho.trim.grav.grd');
system('grdsample -Rgravmodel/thermal.grav.grd gravmodel/rock.grav.grd -Ggravmodel/rock.trim.grav.grd');
system('grdsample -Rgravmodel/thermal.grav.grd gravmodel/sed.grav.grd -Ggravmodel/sed.trim.grav.grd');

% combine interfaces
system(['grdmath gravmodel/flex.moho.trim.grav.grd gravmodel/rock.trim.grav.grd ADD gravmodel/thermal.grav.grd ADD gravmodel/sed.trim.grav.grd ADD = SYNTH.grav.grd']);

% read back in
[XgMod,YgMod,ZgMod]=grdread2('SYNTH.grav.grd');
