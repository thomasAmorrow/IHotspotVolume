%% Hotspot Volume 

% This script calculates along-track volumetric contributions from a
% hotspot trace using bathmetry, age data, and gravity model inputs.

% Bathymetry bounds should provide an adequate buffer around the hotspot
% track such that flexure calculations and gravity models (FFT: subject to
% ringing noise) do not influence the near-hotspot result. Suggested bounds
% are at least 10 degrees long/lat from the edges of the hotspot expression.

% Age data should be placed in AGES.txt following the format in the header
% file.

% Remaining inputs are clipped from global grids using the boundaries of
% the chosen topo file.

% Developed with GMT version 6.1.0_58a5b88_2019.09.08, NOT using GMT-MATLAB API, 
% using OLD gmt script format (GMT 5.x.x compatible)

% For this script to work, calling gmt functions from MATLAB must work.
% 		As a simple test, if the command 
% 			system(['grdsample grdfile -Rgrdfile -I2m+e -G grdfile.trim']);
% 		successfully executes the grdsample function from GMT and returns a 
%		usable grdfile, this script should work

%% ----- REVISION HISTORY -----
% First draft - TMorrow - Aug 03 2018
% Orthogonal interp added - TMorrow Sept 08 2019
% Separating out functions and adding iterative capacity - TMorrow 25 Sep 2019
% Code cleanup - TMorrow 14 Sept. 2020
%  ----- ---------------- -----

%% Set dependencies and libraries

    % clear workspace
    clear
	close all

    % set library paths
	addpath('./dependencies/');
    addpath('./');
    setenv LD_LIBRARY_PATH usr/lib/x86_64-linux-gnu/

%%%% CODE %%%%

%% -- Inputs and paths

    % grid file
    grdname = ' ';
    grdfile=[grdname '.grd'];
    
    % AGES file
	AGES=dlmread('AGES.txt');
    
    % seafloor age file
    SFLagegrd='../../global/infl.age.3.6.grd';
    
    % sediment thickness file
    SEDthckgrd='../../global/sedthick_world_v2.grd';
    
    % WGM FAA file
    WGMFAAgrd='../../global/WGM2012_Freeair_ponc_2min_360.grd';
    
    % adjust AGES data - uncomment if needed
    %AGES(AGES(:,1)<0,1)=AGES(AGES(:,1)<0,1)+360;
    
    % northeast bathy corner
    %   north             east
    NE=[max(AGES(:,2))+10 max(AGES(:,1))+10];
    % southwest bathy corner
    %   south             west
    SW=[min(AGES(:,2))-10 min(AGES(:,1))-10];

    % filter parameters for RR separation
    minW  = 800;  % Minimum filter width candidate for ORS  (e.g., 60) in km
    maxW  = 800; % Maximum filter width candidate for ORS  (e.g., 600) in km
    intW  = 50;  % Filter width step (e.g., 20) in km
    level = 300; % step for base contour calculations
    subaq = 1; % set flag to 1 if all hotspot is underwater
    mask = 1; % set flag to 1 if there are prominent regions that need masking on the map

	% Flexure model inputs
	rho_c =      2800;  % crust density (kg/m3)
	rho_w =      1035;  % water density (kg/m3)
	rho_m =      3300;  % mantle density (kg/m3)
	rho_i =      2400;  % infill density (kg/m3)
    rho_u =      3000;  % underplating density (kg/m3)
	E =          2e23;  % Young's modulus
	v =          0.25;  % Poisson's ratio
	g =          9.8 ;  % Gravity accel (m/s2)
	cr_thck =    7000;  % crust thickness (m)
	T_elas =     400 ;  % elastic isotherm (C)
	T_mant =     1300;  % mantle temperature (C)
	kappa =      1e-6;  % thermal diffusivity

% ---

%% -- (1) Get grid file

    % create CORNERS file for WGET_BATHY
	system('rm CORNERS.xy');
	fileID=fopen('CORNERS.xy','w');
	formatspec='%s \n';
	cornerline=[grdname ',' num2str(NE(2)) ',' num2str(NE(1)) ',' num2str(SW(2)) ',' num2str(SW(1))];
	fprintf(fileID,formatspec,cornerline);
	fclose(fileID);

    % retrieve WGET Script
    system('cp ../../dependencies/WGET_BATHY.sh ./');
    
	% get grid
	system('sh WGET_BATHY.sh');
    
    % resample high res to 2m
    system(['grdsample ' grdfile ' -R' grdfile ' -I2m+e -G' grdfile ]);

    % read in final grid
	[X,Y,Z]=grdread2(grdfile);

% ---

%% -- (2) Mask interfering regions of the map for RR-sep and gravity

	if mask==1
		[XpolyM,YpolyM,INPmask]=IHotVol_PickMask(grdfile);
	end
    
% ---

%% -- (3) Optimized Residual Separation

    % local copy of RR-Sep.sh
    if mask==1
        system('cp ../../dependencies/RR-Sep_mask.sh ./');
        
        system('cp ../../dependencies/RR-Sep-single_mask.sh ./');
    else
        system('cp ../../dependencies/RR-Sep.sh ./');
       
        system('cp ../../dependencies/RR-Sep-single.sh ./');
    end

	% run ORS
	[ORS_L,region]=IHotVol_ORS(grdfile,X,Y,Z,minW,maxW,intW,level,mask);   
    ORS=dlmread('ORStable.txt');
    disp(['Regional/residual separation complete! ORS optimal filter wavelength ' num2str(ORS_L(1)) 'km']);

% ---

%% -- (4) Generate hotspot age track and clipped age sub-grids

	[HSPT_TRK,AGES,pA,TMPTRK]=IHotVol_Track(AGES,X,Y,Z,grdfile,SFLagegrd);

% ---

%% -- (5) Generate initial topographic load grid

	[Xpoly,Ypoly,INP]=IHotVol_PickEdifice(grdfile,AGES);
    
% ---

%% -- (6) Initial flexure solve (Airy)

	IHotVol_Airy(grdfile,rho_c,rho_w,rho_m)

% ---

%% -- (7) Flexure calculation
	
	ii=1;
	[Xflx,Yflx,Zflx]=IHotVol_Flexure([grdfile '_edifice.grd'],rho_c,rho_w,rho_m,rho_i,T_elas,T_mant,kappa,HSPT_TRK,grdfile,ii);

% ---

%% -- (8) Gravity forward model

	% set paths and filenames
	mkdir gravmodel
	sedcutgrdfile='gravmodel/sedcut.DENAN.plusone.grd';
	subaerialgrdfile='gravmodel/subair.DENAN.grd';
	mohoflexgrdfile=['flexure.DENAN.' num2str(ii) '.grd'];
	denangrdfile=['gravmodel/' grdfile '.DENAN.grd'];
	edificegrdfile=[grdfile '_edifice.grd'];

	% make zerogrd for subbing NaNs
	system(['grdmath ' grdfile ' 0 MUL = zerogrdfile.grd']);
	system(['grdsample ' grdfile ' -R' edificegrdfile ' -G' denangrdfile]);

	% DENAN bathy
	system(['grdmath ' denangrdfile ' 0 DENAN = ' denangrdfile ]);

	% Subaerial part of grid
	system(['grdmath ' grdfile ' 0 GT ' grdfile ' MUL = gravmodel/subair.grd']); % TODO look at this
	system(['grdmath gravmodel/subair.grd 0 DENAN = ' subaerialgrdfile ]);

	% trim and DENAN sediment grid for gravity
	system(['grdsample ' SEDthckgrd ' -R' grdfile ' -Ggravmodel/sedcut.' grdfile]);
	system(['grdmath gravmodel/sedcut.' grdfile ' 0 DENAN 1 ' grdfile ' ADD  = ' sedcutgrdfile ]);
	system(['grdmath ' sedcutgrdfile ' 0 DENAN = ' sedcutgrdfile ]);

	% generate forward gravity model
	[XgMod,YgMod,ZgMod]=IHotVol_GravForward(grdfile,denangrdfile,edificegrdfile,sedcutgrdfile,subaerialgrdfile,mohoflexgrdfile,rho_c,rho_w,rho_m,rho_i,kappa,INP,ORS_L,subaq);
		
% ---

%% -- (9) FAA

	% get residual
    [XResG,YResG,ZResG]=IHotVol_FAAgetResidual(ORS_L,WGMFAAgrd,mask,subaq);
    
    % polygon select for residual determination
    INPold=INP;
    [XResGmesh,YResGmesh]=meshgrid(XResG,YResG);
    INP=inpolygon(XResGmesh,YResGmesh,Xpoly,Ypoly);
    
    % gravity residual along HSPT_TRK
    grav_resid=sqrt(double(sum(sum((INP.*ZResG).^2))));
    
% ---

%% -- (10) Underplating

    [Xg,Yg,finaltopoinverse]=IHotVol_Underplating(Xflx,Yflx,Zflx,XResG,YResG,ZResG,ii,grdfile,ORS_L,1e-5);
    
% ---

%% -- (11) Iterative flexure/underplating revisions

	% convergence criterion setup
	resisDiff=[];
	residDiff(1)=1e10;

	% new load from compensated topography
	while residDiff(ii)>0.0001
		
		ii=ii+1

		% resample/fit
		system(['grdsample Uplate.' num2str(ii-1) '.grd -R' grdfile '_edifice.grd -GUplate.' num2str(ii-1) '.grd']);
		system(['grdmath Uplate.' num2str(ii-1) '.grd  Uplate.' num2str(ii-1) '.grd LOWER SUB = Uplate.' num2str(ii-1) '.grd']);
		
		% load reduction correction
		system(['grdmath ' grdfile '_edifice.grd ' num2str(rho_c-rho_w) ' MUL Uplate.' num2str(ii-1) '.grd 1000 MUL ' num2str(rho_u-rho_m) ' MUL ADD ' num2str(rho_c-rho_w) ' DIV 0 DENAN = ' grdfile '_edifice.' num2str(ii) '.grd']);
		
		% new flexure profile
		[Xflx,Yflx,Zflx]=IHotVol_Flexure([grdfile '_edifice.' num2str(ii) '.grd'],rho_c,rho_w,rho_m,rho_i,T_elas,T_mant,kappa,HSPT_TRK,grdfile,ii);
		mohoflexgrdfile=['flexure.DENAN.' num2str(ii) '.grd'];
		
		
		if ii==2
		    system(['grdsample ' grdfile ' -R' grdfile '_edifice.' num2str(ii) '.grd -G' grdfile 'sampleiter.grd']); 
		    system(['grdmath ' grdfile '_edifice.' num2str(ii) '.grd ' grdfile '_edifice.' num2str(ii) '.grd MEAN ' grdfile 'sampleiter.grd MEAN SUB SUB 0 DENAN = ' grdfile '_edifice.' num2str(ii) '.grd']);
		else
		    system(['grdmath ' grdfile '_edifice.' num2str(ii) '.grd ' grdfile '_edifice.' num2str(ii) '.grd MEAN ' grdfile '_edifice.' num2str(ii-1) '.grd MEAN SUB SUB 0 DENAN = ' grdfile '_edifice.' num2str(ii) '.grd']);
		end

		% forward gravity model		
		[XgMod,YgMod,ZgMod]=IHotVol_GravForward(grdfile,denangrdfile,edificegrdfile,sedcutgrdfile,subaerialgrdfile,mohoflexgrdfile,rho_c,rho_w,rho_m,rho_i,kappa,INP,ORS_L,subaq);

		% new residual
		[XResG,YResG,ZResGNEW]=IHotVol_FAAgetResidual(ORS_L,WGMFAAgrd,mask,subaq);
		
		% gravity residual along HSPT_TRK
		grav_residNEW=sqrt(double(sum(sum((INP.*ZResGNEW).^2))));
		
		% assess convergence
		residDiff(ii)=abs(grav_resid(ii-1)-grav_residNEW);
		grav_resid(ii)=grav_residNEW;
		disp(['current iteration residual change:' num2str(residDiff(ii))])
		
		% calculate underplating
		[Xg,Yg,finaltopoinverse]=IHotVol_Underplating(Xflx,Yflx,Zflx,XResG,YResG,ZResG,ii,grdfile,ORS_L,1e-5);
		
	end


	% final flexure calculation using only the compensated edifice
	system(['grdsample ' grdfile '_edifice.grd -R' grdfile '_edifice.' num2str(ii) '.grd -G' grdfile '_edifice.flexsample.grd']);
	system(['grdmath ' grdfile '_edifice.flexsample.grd 0 GT ' grdfile '_edifice.' num2str(ii) '.grd MUL 0 DENAN = ' grdfile '_edifice.FINAL.grd']);
	[Xflx,Yflx,Zflx]=IHotVol_Flexure([grdfile '_edifice.' num2str(ii) '.grd'],rho_c,rho_w,rho_m,rho_i,T_elas,T_mant,kappa,HSPT_TRK,grdfile,ii);

% ---

%% -- (12) Sample final output grids 

	IHotVol_SampleVolGrids(grdfile,ii,HSPT_TRK,mask);

% ---

%% -- (13) Import cross sections

	% this is best done by hand at this point.
	% import the text file as a matrix 
	% it should default to a variable named 'grdname'

	% --**-- IMPORTANT --**--
	% in Calculate Volumes below, the variable passed to the 
	% VolumeSlices function must match the name of the imported cross-sections
	% --**-- --------- --**--

% ---

%% -- (14) Calculate volumes

	close all
	VOL=[];
	
	% volume calculation
	[VOL,Crosses]=IHotVol_VolumeSlices(grdname,HSPT_TRK);

	% generate age-volume plot
	VOL=IHotVol_AgeVolPlot(VOL,Crosses,pA);

% ---

%% -- (15) Spectral analysis

	PEAKS=IHotVol_Spectral(VOL,grdfile);

% ---

%% -- (16) Save workspace

	save completed
	
% ---
