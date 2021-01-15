function [Xpoly,Ypoly,INP]=IHotVol_PickEdifice(grdfile,AGES)

%% Hand pick edifice from residual-separated bathymetry

% load residual grid and contour
[Xres,Yres,Zres]=grdread2([grdfile '_residual.grd']);

% zero out negative anomalies
Zres(Zres<=0)=0;

close all

% PICKING ------------ %

% contour residual bathy and plot age data
figure
C=contourf(Xres,Yres,Zres,20);
title('Residual bathymetry');
axis equal
hold on
plot(AGES(:,1),AGES(:,2),'ro-');
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
system(['grd2xyz ED.grd | xyz2grd -G' grdfile '_edifice.grd -R' grdfile '.trim']);

% de-nan the edifice grid (0 outside)
system(['grdmath ' grdfile '_edifice.grd 0 DENAN = ' grdfile '_edifice.grd']);

% write out the in-polygon grid
[Xg,Yg,Zg]=grdread2([grdfile '_edifice.grd']);
grdwrite2(Xg,Yg,double(INP),'INP.grd');

% contour edifice grid
figure
C=contourf(Xg,Yg,Zg,20);
title('Edifice');
axis equal;	

% final display
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]);
disp('You picked this material as your volcanic edifice');
drawnow;
