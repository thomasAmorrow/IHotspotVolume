function [HSPT_TRK,AGES,pA,TMPTRK]=IHotVol_Track(AGES,X,Y,Z,grdfile,SFLagegrd)

%% Upload hotspot age track and determine least squares polynomial fit for both age and location

% sort age data
AGES=sortrows(AGES,1);

% Lat/long relocation ------------ %
x_track=[];
y_track=[];
age_track=[];
Xtrk=AGES(:,1);
Ytrk=AGES(:,2);
Lage=AGES(:,3);

% set polynomial maximum order
if length(AGES(:,1))<6
    pmax=length(AGES(:,1))-1;
else
    pmax=6;
end

close all

% create inspection map
contourf(X,Y,Z,20);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]);
axis equal
hold on
contour(X,Y,Z,20);
plot(AGES(:,1),AGES(:,2),'ko');

% generate polynomial fits
for n=1:pmax
    [p S]=polyfit(Xtrk,Ytrk,n);
    plot(min(AGES(:,1)):0.1:max(AGES(:,1)),polyval(p,min(AGES(:,1)):0.1:max(AGES(:,1))),'Linewidth',2,'DisplayName',num2str(n));
    Fit(n)=S.normr;
    hold on
end

% visual inspection for polynomials
plot(AGES(:,1),AGES(:,2),'ko');
legend('show');
drawnow
disp('What order polynomial best follows the hotspot?');
numFit=input('>:');

% polynomial fitting
[p S]=polyfit(Xtrk,Ytrk,numFit);

% orthogonally minimized (least squares) relocation of age data
for ii=1:length(Xtrk)
	[x_track(ii),y_track(ii)]=closepoint(AGES(ii,1),AGES(ii,2),p);
end
TMPTRK=[x_track' y_track' AGES(:,3)];
TMPTRK=sortrows(TMPTRK,1);

%*** UNCOMMENT THIS LINE FOR HOTSPOTS CROSSING +/-180 longitude ***%
%TMPTRK(TMPTRK(:,1)<0,1)=TMPTRK(TMPTRK(:,1)<0,1)+360;
%*** UNCOMMENT THIS LINE FOR HOTSPOTS CROSSING +/-180 longitude ***%

% output to text file
writematrix(TMPTRK,'track.txt');

% generate track with GMT
system(['cat track.txt | awk ''{print $1,$2}'' | grdtrack -G' grdfile ' -C2k/1/20 -Ar > track_cross.txt']);
cmdstr='cat track_cross.txt | awk ''{if ($3==0) print $1,$2}'' > track20k.txt';
system(cmdstr);

% read in to MATLAB
TRK=dlmread('track20k.txt');

%*** UNCOMMENT THIS LINE FOR HOTSPOTS CROSSING +/-180 longitude ***%
%TRK(TRK(:,1)<0,1)=TRK(TRK(:,1)<0,1)+360;
%*** UNCOMMENT THIS LINE FOR HOTSPOTS CROSSING +/-180 longitude ***%

TMPTRK=sortrows(TMPTRK,3);


% get distance-along for age interpretation
DIST=sqrt((TRK(:,1)-TMPTRK(1,1)).^2+(TRK(:,2)-TMPTRK(1,2)).^2);

if find(DIST==min(DIST))<0.5*length(TMPTRK(:,1))
    TRK(:,3)=20.*(0:length(TRK(:,1))-1)';
else
    TRK=flipud(TRK);
    TRK(:,3)=20.*(0:length(TRK(:,1))-1)';
end

% distance fit based on same order polynomial from lat-long fit
pDA=polyfit(TRK(:,1),TRK(:,3),numFit);
TMPTRK(:,4)=polyval(pDA,TMPTRK(:,1));

close all

% plot age-distance
plot(TMPTRK(:,4),TMPTRK(:,3),'rs')

% generate polynomial fits
for n=1:5
    [pA S]=polyfit(TMPTRK(:,4),TMPTRK(:,3),n);
    plot(min(TMPTRK(:,4)):1:max(TMPTRK(:,4)),polyval(pA,min(TMPTRK(:,4)):1:max(TMPTRK(:,4))),'Linewidth',2,'DisplayName',num2str(n));
    Fit(n)=S.normr;
    hold on
end

% visual inspection for polynomials
plot(TMPTRK(:,4),TMPTRK(:,3),'rs')
legend('show');
drawnow
disp('What order polynomial best follows the age/distance relationship?');
numFit=input('>:');

% polynomial fitting
[pA S]=polyfit(TMPTRK(:,4),TMPTRK(:,3),numFit);

TRK(:,4)=polyval(pA,TRK(:,3));

% build final output
HSPT_X=[]; 
HSPT_Y=[];
HSPT_A=[];
HSPT_X=TRK(:,1);
HSPT_Y=TRK(:,2);
HSPT_A=TRK(:,4);
HSPT_L=TRK(:,3);

% concatenate the results into a single array
HSPT_TRK=[HSPT_X HSPT_Y HSPT_A];

% directory for age grid subsamples
!mkdir cutagetmp

% Read age grid file, cut to bathy bounds
system(['grdsample ' SFLagegrd ' -Gcutagetmp/age.' grdfile ' -R' grdfile]);

% get weighted-average seafloor age
for ii=1:length(HSPT_TRK(:,1))
	clear tmpX tmpY tmpZ WTG
	% 0.5-degree circle cut from age grid
	system(['grdcut cutagetmp/age.' grdfile ' -Sn' num2str(HSPT_TRK(ii,1)) '/' ...
		num2str(HSPT_TRK(ii,2)) '/5d -Gcutagetmp/cut.age.' num2str(ii) '.grd']);
    
    % read in local age grid
    [tmpX tmpY tmpZ]=grdread2(['cutagetmp/cut.age.' num2str(ii) '.grd']);

	% weighting
	for jj=1:length(tmpX(1,:));
		for kk=1:length(tmpY(1,:));
			WTG(kk,jj)=1/distance(HSPT_TRK(ii,2),HSPT_TRK(ii,1),tmpY(kk),tmpX(jj));
		end
	end

	% append to hotspot age data
	HSPT_TRK(ii,4)=nansum(reshape(WTG,numel(WTG),1).*reshape(tmpZ,numel(tmpZ),1))/...
	nansum(reshape(WTG,numel(WTG),1));

	clear tmpX tmpY tmpZ WTG

end

% append loading age
HSPT_TRK=[HSPT_X HSPT_Y HSPT_A HSPT_TRK(:,4) HSPT_L];
HSPT_TRK(:,3)=HSPT_TRK(:,3)-min(HSPT_TRK(:,3));

% make sure no duplicates
HSPT_TRK=unique(HSPT_TRK,'rows');
% ---
