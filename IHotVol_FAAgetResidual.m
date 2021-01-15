function[XResG,YResG,ZResG]=IHotVol_FAAgetResidual(ORS_L,WGMFAAgrd,mask,subaq)

% determines FAA anomaly and residual from synthetic gravity approximation

% Sample WGM data
system(['grdsample ' WGMFAAgrd ...
    	' -RSYNTH.grav.grd -Ggravmodel/faa.cut.grd']);

% subtract forward model to get residual
system('grdmath gravmodel/faa.cut.grd SYNTH.grav.grd SUB = gravmodel/UPLATE.grav.grd');

% remove mean
system(['grdmath gravmodel/faa.cut.grd gravmodel/faa.cut.grd MEAN SUB = gravmodel/faa.cut.grd']);

% remove mean of forward model
system(['grdmath SYNTH.grav.grd SYNTH.grav.grd MEAN SUB = SYNTH.grav.grd']);

% low-pass filter underplating grid at 1/2 edifice wavelength
system(['grdfilter gravmodel/UPLATE.grav.grd -Fc20k -D1 -Ggravmodel/UPLATE.grav.grd']);
!cp gravmodel/UPLATE.grav.grd gravmodel/UPLATE.lowpass.grav.grd

% if masked, also mask result
if mask==1
    system(['grdsample MASK.grd -Rgravmodel/UPLATE.grav.grd -GMASK.grav.grd']);
    system(['cp gravmodel/UPLATE.grav.grd gravmodel/backup.UPLATE.grav.grd']);
    system(['grdmath MASK.grav.grd gravmodel/UPLATE.grav.grd MUL = gravmodel/UPLATE.grav.grd']);
end

% read in results
[XResG,YResG,ZResG]=grdread2('gravmodel/UPLATE.grav.grd');

% pad edges to help with inversion
ZResG(:,1)=ZResG(:,2);
ZResG(:,end)=ZResG(:,end-1);
ZResG(end,:)=ZResG(end-1,:);
ZResG(1,:)=ZResG(2,:);
