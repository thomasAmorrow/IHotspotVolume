function[]=IHotVol_SampleVolGrids(grdfile,ii,HSPT_TRK,mask)
    
    
% formatting grids 
system(['grd2xyz Uplate.' num2str(ii) '.grd | xyz2grd -GUplate.' num2str(ii) '.grd -R' grdfile]);
system(['grd2xyz flexure.DENAN.' num2str(ii) '.grd | xyz2grd -Gflexure.DENAN.' num2str(ii) '.grd -R' grdfile]);
system(['grd2xyz ' grdfile '_edifice.grd | xyz2grd -G' grdfile '_edifice.grd -R' grdfile]);
system(['grd2xyz INP.grd | xyz2grd -GINP.grd -R' grdfile]);

% sample the inpolygon grid again
system(['grdsample INP.grd -RUplate.' num2str(ii) '.grd -GINP.sampled.grd']);

% determine underplating within the polygon region
system(['grdmath INP.sampled.grd Uplate.' num2str(ii) '.grd MUL = Uplate.INP.' num2str(ii) '.grd']);

%% -- Secondary masking (if needed)
quickcontgrd(['Uplate.' num2str(ii) '.grd'],20);
quickcontgrd(['flexure.DENAN.' num2str(ii) '.grd'],20);

disp('Do the results need another masking? 1:y 0:n')
mask2=input('>:');

% additional masking
if mask2==1
    [XpolyM,YpolyM,INPmask]=IHotVol_PickMask2(['Uplate.' num2str(ii) '.grd']);

    
    system(['grdmath flexure.DENAN.' num2str(ii) '.grd MASK2.grd MUL 0 NAN = flexure.DENAN.' num2str(ii) '.grd']);
    system(['grdmath Uplate.' num2str(ii) '.grd MASK2.grd MUL 0 NAN = Uplate.' num2str(ii) '.grd']);


end

% ---

% remove conflicting files
system(['rm ' grdfile 'Xes.txt']);

% sample resultant grids
system(['cat track.txt | awk ''{print $1,$2}'' | grdtrack -R' grdfile ' -V -G' ...
    grdfile '_edifice.grd -Gflexure.DENAN.' num2str(ii) '.grd -GUplate.' num2str(ii) '.grd -GUplate.INP.' num2str(ii) '.grd -C1000k/4k/4k -Ar > ' grdfile 'Xes.txt']);

% ----
