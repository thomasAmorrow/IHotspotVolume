function [Xpoly,Ypoly,INP]=IHotVol_PickEdifice(grdfile,AGES)

%% Hand pick edifice from residual-separated bathymetry

% load residual grid and contour
[Xres,Yres,Zres]=grdread2([grdfile '_residual.grd']);

close all

% PICKING ------------ %

% contour residual bathy
figure
C=contourf(Xres,Yres,Zres,12);
title('Residual bathymetry');
axis equal
hold on
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]);
pause(5);

disp(['I cleaned this up as best I can. Please pick a polygon around ' ...
'the edifice of interest']);

% select edifice bounds - will zero residual outside polygon
[Xpoly,Ypoly]=ginput;

% ------------ %

close all

% meshing the X and Y values for polygon 
[Xg,Yg]=meshgrid(Xres,Yres);

% find grid values within polygon
INP=inpolygon(Xg,Yg,Xpoly,Ypoly);
Zed=Zres.*INP;
Edifice=Zed;
Zed(Zed==0)=nan;

% output edifice grid
grdwrite2(Xres,Yres,Edifice,'ED.grd')
grdwrite2(Xres,Yres,double(INP),'INP.grd')

% reformat grid (output from grdwrite2 is deprecated and does not cooperate with grdflexure)
system(['grd2xyz ED.grd | xyz2grd -G' grdfile '_edifice.grd -R' grdfile]);

% de-nan the edifice grid (0 outside)
system(['grdmath ' grdfile '_edifice.grd 0 DENAN = ' grdfile '_edifice.grd']);

% read back in
[Xg,Yg,Zg]=grdread2([grdfile '_edifice.grd']);

% contour edifice grid
figure
C=contourf(Xg,Yg,Zg,20);
title('Edifice');
axis equal;	

% show output
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]);
disp('You picked this material as your volcanic edifice');
drawnow;
