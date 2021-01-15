function [ORS_L,region] = IHotVol_ORS(grdfile,X,Y,Z,minW,maxW,intW,level,mask)

%% Optimized residual separation (Wessel, 2016)

% sub-region of map
region  = ['-R' num2str(ceil(min(X))+0.05) '/' num2str(floor(max(X))-0.05) ... 
			'/' num2str(ceil(min(Y))+0.05) '/' num2str(floor(max(Y))-0.05)];
			
% mask requires a different RR-Sep script			
if mask==1        
    system(['cp ' grdfile ' backupgrdfile_ORS.grd']);
    system(['grdmath MASK.grd ' grdfile ' MUL 0 DENAN = ' grdfile]);
    cmdStr = ['sh ./RR-Sep_mask.sh ' grdfile ' ' region ' ' num2str(minW) ...
        ' ' num2str(maxW) ' ' num2str(intW) ' ' num2str(level)];
    system(cmdStr);
else
    cmdStr = ['sh ./RR-Sep.sh ' grdfile ' ' region ' ' num2str(minW) ...
        ' ' num2str(maxW) ' ' num2str(intW) ' ' num2str(level)];
    system(cmdStr);
end

% read ORS table (EXISTING)
ORS=dlmread('ORStable.txt');

% get ORS max
ORS_L=ORS(ORS(:,5)==max(ORS(:,5)),:);

% 0. System defaults (Change if you know what you are doing):
dim_dist=2;		% How we compute distances on the grid [Flat Earth approximation]
dim_sectors=8;		% Number of sectors to use [8]
dim_filter='m';		% Primary filter [m for median]
dim_quantity='l';		% Secondary filter [l for lower]
dim_smooth_type='m';	% Smoothing filter type [m for median]
dim_smooth_width=50;	% Smoothing filter width, in km [50]
orsout='ORSout';

if mask==1
    system(['cp backupgrdfile_ORS.grd ' grdfile ])
end

% 0. System defaults (Change if you know what you are doing):
dim_dist=2;		% How we compute distances on the grid [Flat Earth approximation]
dim_sectors=8;		% Number of sectors to use [8]
dim_filter='m';		% Primary filter [m for median]
dim_quantity='l';		% Secondary filter [l for lower]
dim_smooth_type='m';	% Smoothing filter type [m for median]
dim_smooth_width=50;	% Smoothing filter width, in km [50]
orsout='ORSout';


% final filtering
disp('Filter wavelength found, performing final separation')

if mask==1
    cmdStr = ['sh ./RR-Sep-single_mask.sh ' grdfile ' ' region ' ' num2str(ORS_L(1)) ...
        ' ' num2str(ORS_L(1)) ' ' num2str(intW) ' ' num2str(level)];
    system(cmdStr);
else
    cmdStr = ['sh ./RR-Sep-single.sh ' grdfile ' ' region ' ' num2str(ORS_L(1)) ...
        ' ' num2str(ORS_L(1)) ' ' num2str(intW) ' ' num2str(level)];
    system(cmdStr);
end

% subtract ORS regional from observed
system(['cp ./finalRR/resid.grd ' grdfile '_residual.grd']);
system(['grdsample ' grdfile ' -R' grdfile '_residual.grd -Gtmpregion.grd']);
system(['grdmath tmpregion.grd ' grdfile '_residual.grd SUB = ' grdfile '_regional.grd']);

% ----
