function []=IHotVol_Airy(grdfile,rho_c,rho_w,rho_m)

%% Airy isostatic underplating estimate

% Airy compensation calculated directly from edifice grid
system(['grdmath ' num2str(rho_c) ' ' num2str(rho_w) ' SUB ' num2str(rho_m) ' ' ...
        num2str(rho_c) ' SUB DIV ' grdfile '_edifice.grd MUL NEG = ' grdfile ...
		'_airy_compensation.grd']);

