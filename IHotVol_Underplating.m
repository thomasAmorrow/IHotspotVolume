function[Xg,Yg,finaltopoinverse,rms]=IHotVol_Underplating(Xfl,Yfl,Zfl,Xg,Yg,Zg,ii,grdfile,ORS_L,criterio)

%% -- Invert gravity for underplating

	% Modified from 3DINVER.M: A MATLAB program to invert the gravity anomaly over a 3-D horizontal density interface by Parker-Oldenburg's algorithm
	% David Gomez Ortiz and Bhrigu N P Agarwal

	% Iterative process of Oldenburg(1974)
	% using the Parker's (1973) equation
	% Given a gravity map as input (bouin)in mGal
	% as well as density contrast (contrast) in gr/cm3
	% and the mean source depth (z0) in Km
	% the program computes the topography of the interface (topoout) in Km 

	% density contrast is positive (e. g. mantle-crust)
	% mean depth reference is positive downwards(z0)
	% It is also necessary to include map length in Km (longx and longy)
	% considering a square area, the number of rows or columns (numrows and numcolumns)
	% and the convergence criterion (criterio) in Km

	% The maximum number of iterations is 10

	% The program computes also the gravity map (bouinv) due to the computed relief(forward modelling)
	% The cut-off frequencies (SH and WH) for the filter in order to achieve
	% convergence have to be also previously defined (as 1/wavelength in Km)

%% code

    % Mirror Z grid, X, Y to get larger grid, resample down to square
    Zg1=flipud(Zg);
    Zg2=fliplr(Zg);
    Zg3=rot90(Zg,2);
    Zstitch=[Zg3 Zg1 Zg3; Zg2 Zg Zg2; Zg3 Zg1 Zg3];
    Xstitch=[Xg-abs(Xg(1)-Xg(end))-abs(Xg(1)-Xg(2)) Xg Xg+abs(Xg(1)-Xg(end))+abs(Xg(end-1)-Xg(end))];
    Ystitch=[Yg-abs(Yg(1)-Yg(end))-abs(Yg(1)-Yg(2)) Yg Yg+abs(Yg(1)-Yg(end))+abs(Yg(end-1)-Yg(end))];

	% box is next power of 2 larger than largest dimension of input
    boxEl=max([2^nextpow2(numel(Xg)) 2^nextpow2(numel(Yg))]);
    
    if boxEl==2^nextpow2(numel(Xg)) && boxEl>3*numel(Yg)
        boxEl=max([2^(nextpow2(numel(Xg))-1)+0.5*(2^(nextpow2(numel(Xg)))-2^(nextpow2(numel(Xg))-1)) 2^(nextpow2(numel(Yg))-1)+0.5*(2^(nextpow2(numel(Yg)))-2^(nextpow2(numel(Yg))-1))]);
    end
    if boxEl==2^nextpow2(numel(Yg)) && boxEl>3*numel(Xg)
        boxEl=max([2^(nextpow2(numel(Xg))-1)+0.5*(2^(nextpow2(numel(Xg)))-2^(nextpow2(numel(Xg))-1)) 2^(nextpow2(numel(Yg))-1)+0.5*(2^(nextpow2(numel(Yg)))-2^(nextpow2(numel(Yg))-1))]);
    end
    
    Xind=round(numel(Xg)/2)+numel(Xg);
    Yind=round(numel(Yg)/2)+numel(Yg);
    
    Xsq=Xstitch((Xind-boxEl/2):(Xind+boxEl/2)-1);
    Ysq=Ystitch((Yind-boxEl/2):(Yind+boxEl/2)-1);
    
    Zsq=Zstitch((Yind-boxEl/2):(Yind+boxEl/2-1),(Xind-boxEl/2):(Xind+boxEl/2-1));

    % square sampled, power of 2 grid is ready now
    Xg=Xsq;
    Yg=Ysq;
    Zg=Zsq;
	Zg=double(Zg);

	% size of box
	numrows=numel(Yg);              
	numcolumns=numel(Xg);           
   
   	% dimension box 
	longx=deg2km(distance(Yg(numel(Yg)/2),Xg(1),Yg(numel(Yg)/2),Xg(end)));              %e.g. 666 Km data length in X direction
	longy=deg2km(distance(Yg(1),Xg(numel(Xg)/2),Yg(end),Xg(numel(Xg)/2)));            %e.g. 666 Km data length in Y direction
	contrast=0.3;          %density contrast value in g/cm3, e. g. 0.4
	z0=11+abs(nanmin(nanmin(Zfl)))/1e3;               %mean reference depth in Km, e. g. 30 Km
	WH=1/100;                %smaller cut-off frequency (1/wavelength in Km) e. g. 0.01, i.e. 1/100 Km
	SH=1/80; %(2*ORS_L(1)-0.1*ORS_L(1));                %greater cut-off frequency (1/wavelength in Km) e. g. 0.012, i.e. 1/83.3 Km

	truncation=0.1;          %truncation window data length in % per one, i.e, 10%

	%input gravity map is demeaned
	bou=Zg;
    meangravity=nanmean(nanmean(Zg));
    bou=bou-meangravity;  
    
	fftbou=fft2(bou);   %fft2 computes the 2-D FFT of a 2-D matrix (bou, in this case)

	%A cosine Tukey window with a truncation of 10% default is applied
	wrows = tukeywin(numel(Yg),truncation);  %this computes a 1-D cosine Tukey window of the same length as the original matrix input rows and with the truncation defined by the variable 'truncation'
	wcolumns = tukeywin(numel(Xg),truncation);  %this computes a 1-D cosine Tukey window of the same length as the original matrix input columns and with the truncation defined by the variable 'truncation'

	w2 =wrows* wcolumns.'; %this generates a 2-D cosine Tukey window multipliying the 1-D windows
	bou=bou.*w2;   %the original gravity input matrix, previously demeaned, is multiplied by the cosine window

	mapabou=bou;   %the original gravity input matrix after demeaning is transposed

	fftbou=fft2(mapabou);  %the 2-D FFT of the gravity input matrix is computed after demeaning
	spectrum=abs(fftbou(1:numrows/2, 1:numcolumns/2));  %this computes the amplitude spectrum

	%the matrix with the frequencies of every harmonic is computed
	for f=1:abs((numrows/2)+1)
	   for g=1:abs((numcolumns/2)+1)
	      frequency(f, g)=sqrt(((f-1)/longx)^2+((g-1)/longy)^2);
       end
    end
    
	%the matrix of the negative frequencies is also computed
	frequency2=fliplr(frequency);
	frequency3=flipud(frequency);
	frequency4=fliplr(flipud(frequency));
	entero=round(numcolumns/2);
	if ((numcolumns/2)- entero)==0
	   frequency2(:,1)=[];
	   frequency3(1,:)=[];
	   frequency4(:,1)=[];
	   frequency4(1,:)=[];
	   frequencytotal=[frequency frequency2;frequency3 frequency4];
	else
	   frequencytotal=[frequency frequency2;frequency3 frequency4];
	end
	frequencytotal(end,:)=[];
	frequencytotal(:,end)=[];

	frequencytotal=frequencytotal.*(2*pi);  %the frequency (1/wavelength) matrix is transformed to wavenumber (2*pi/wavelength) matrix

	%The iterative process starts here
	%The first term of the series, that is constant, is computed and stored in variable 'constant'
	up=-(fftbou.*(exp((z0)*(frequencytotal))));
	down=(2*pi*6.67*contrast);
	constant=up./down;
	
	%The high-cut filter is constructed
	filter=frequencytotal.*0;      %the filter matrix is set to zero
	frequencytotal=frequencytotal./(2*pi);  
	for f=1:numrows
		for g=1:numcolumns
	        if frequencytotal(f,g)<WH
	            filter(f,g)=1;  
	        elseif frequencytotal(f,g)<SH
	            filter(f,g)=0.5.*(1+cos((((2*pi)*frequencytotal(f,g))-(2*pi*WH))/(2*(SH-WH))));
	        else
	            filter(f,g)=0;
            end
	    end
	end
	constant=constant.*filter;  %the filter is applied to the first term of the series

	topoinverse=real(ifft2(constant));  %the real part of the inverse 2-D FFT of the first term provides the first topography approach stored in the variable 'topoinverse'
	frequencytotal=frequencytotal.*(2*pi);

	%It starts the computation of the second term of the series with a maximum of 10 iterations
   
	topo2=(((frequencytotal.^(1))./(prod(1:2))).*(fft2(topoinverse.^2)));   %the function 'prod' computes the factorial of the number inside the brackets

	topo2=topo2.*filter;   %the filter is applied to the second term of the series

	topo2=constant-topo2; %the new topography approach is computed substracting the first term to the second one

	topoinverse2=real(ifft2(topo2));   % the real part of the 2-D inverse FFT provides the new topography approach in space domain

	diference2=topoinverse2-topoinverse;  % this compute the diference between the two conscutive topography approaches
	diference2=diference2.^2;
	rms2=sqrt(sum(sum(diference2))/(2*(numrows*numcolumns)));   %it computes the rms error between the last two iterations
	iter=2;  %the number of the iteration reached is stored in variable 'iter'
	rms=rms2
	finaltopoinverse=topoinverse2;  %the new topography approach is stored in the variable 'finaltpoinverse'
	if rms2>=criterio   %it determines if the rms error is greater than the convergence criterion
	    topo3=(((frequencytotal.^(1))./(prod(1:2))).*(fft2(finaltopoinverse.^2)))+((((frequencytotal.^(2))./(prod(1:3))).*(fft2(finaltopoinverse.^3))));   %the third term of the series is computed using the new topography approach


	topo3=topo3.*filter; %the filter is applied to the third term of the series

	topo3=constant-topo3;  %the new topography approach is computed
    topoinverse3=real(ifft2(topo3));  %and transformed to space domain
    diference3=topoinverse3-topoinverse2;  %the diference between the last two topography approaches is computed
    diference3=diference3.^2;
    rms3=sqrt(sum(sum(diference3))/(2*(numrows*numcolumns)));  %the new rms error is computed
    rms=rms3
    finaltopoinverse=topoinverse3;  %the new topography approach is stored in variable finaltopoinverse
    iter=3;  %it indicates the iteration step reached
	if rms3>=criterio  %it determines if the convergence criterion has been reached or not
    	topo4=(((frequencytotal.^(1))/(prod(1:2))).*(fft2(topoinverse3.^2)))+((((frequencytotal.^(2))/(prod(1:3))).*(fft2(topoinverse3.^3))))+((((frequencytotal.^(3))/(prod(1:4))).*(fft2(topoinverse3.^4))));  %the fourth term of the series is computed, and so on...

	topo4=topo4.*filter;

	topo4=constant-topo4;
    topoinverse4=real(ifft2(topo4));
    diference4=topoinverse4-topoinverse3;
    diference4=diference4.^2;
    rms4=sqrt(sum(sum(diference4))/(2*(numrows*numcolumns)));
    rms=rms4
    finaltopoinverse=topoinverse4;
    iter=4;
    if rms4>=criterio
         topo5=(((frequencytotal.^(1))/(prod(1:2))).*(fft2(topoinverse4.^2)))+((((frequencytotal.^(2))/(prod(1:3))).*(fft2(topoinverse4.^3))))+((((frequencytotal.^(3))/(prod(1:4))).*(fft2(topoinverse4.^4))))+((((frequencytotal.^(4))/(prod(1:5))).*(fft2(topoinverse4.^5))));

	topo5=topo5.*filter;

    		topo5=constant-topo5;
    		topoinverse5=real(ifft2(topo5));
    		diference5=topoinverse5-topoinverse4;
    		diference5=diference5.^2;
      	rms5=sqrt(sum(sum(diference5))/(2*(numrows*numcolumns)));
         rms=rms5
         finaltopoinverse=topoinverse5;
         iter=5;
         if rms5>=criterio
            topo6=(((frequencytotal.^(1))/(prod(1:2))).*(fft2(topoinverse5.^2)))+((((frequencytotal.^(2))/(prod(1:3))).*(fft2(topoinverse5.^3))))+((((frequencytotal.^(3))/(prod(1:4))).*(fft2(topoinverse5.^4))))+((((frequencytotal.^(4))/(prod(1:5))).*(fft2(topoinverse5.^5))))+((((frequencytotal.^(5))/(prod(1:6))).*(fft2(topoinverse5.^6))));

topo6=topo6.*filter;

    			topo6=constant-topo6;
    			topoinverse6=real(ifft2(topo6));
    			diference6=topoinverse6-topoinverse5;
    			diference6=diference6.^2;
      		rms6=sqrt(sum(sum(diference6))/(2*(numrows*numcolumns)));
            rms=rms6
            finaltopoinverse=topoinverse6;
            iter=6;
            if rms6>=criterio
               topo7=(((frequencytotal.^(1))/(prod(1:2))).*(fft2(topoinverse6.^2)))+((((frequencytotal.^(2))/(prod(1:3))).*(fft2(topoinverse6.^3))))+((((frequencytotal.^(3))/(prod(1:4))).*(fft2(topoinverse6.^4))))+((((frequencytotal.^(4))/(prod(1:5))).*(fft2(topoinverse6.^5))))+((((frequencytotal.^(5))/(prod(1:6))).*(fft2(topoinverse6.^6))))+((((frequencytotal.^(6))/(prod(1:7))).*(fft2(topoinverse6.^7))));

topo7=topo7.*filter;

    		topo7=constant-topo7;
    		topoinverse7=real(ifft2(topo7));
    		diference7=topoinverse7-topoinverse6;
    				diference7=diference7.^2;
      			rms7=sqrt(sum(sum(diference7))/(2*(numrows*numcolumns)));
               rms=rms7
               finaltopoinverse=topoinverse7;
               iter=7;
               if rms7>=criterio
                  topo8=(((frequencytotal.^(1))/(prod(1:2))).*(fft2(topoinverse7.^2)))+((((frequencytotal.^(2))/(prod(1:3))).*(fft2(topoinverse7.^3))))+((((frequencytotal.^(3))/(prod(1:4))).*(fft2(topoinverse7.^4))))+((((frequencytotal.^(4))/(prod(1:5))).*(fft2(topoinverse7.^5))))+((((frequencytotal.^(5))/(prod(1:6))).*(fft2(topoinverse7.^6))))+((((frequencytotal.^(6))/(prod(1:7))).*(fft2(topoinverse7.^7))))+((((frequencytotal.^(7))/(prod(1:8))).*(fft2(topoinverse7.^8))));

topo8=topo8.*filter;
    					topo8=constant-topo8;
    					topoinverse8=real(ifft2(topo8));
    					diference8=topoinverse8-topoinverse7;
    					diference8=diference8.^2;
      				rms8=sqrt(sum(sum(diference8))/(2*(numrows*numcolumns)));
                  rms=rms8
                  finaltopoinverse=topoinverse8;
                  iter=8;
                  if rms8>=criterio
                     topo9=(((frequencytotal.^(1))/(prod(1:2))).*(fft2(topoinverse8.^2)))+((((frequencytotal.^(2))/(prod(1:3))).*(fft2(topoinverse8.^3))))+((((frequencytotal.^(3))/(prod(1:4))).*(fft2(topoinverse8.^4))))+((((frequencytotal.^(4))/(prod(1:5))).*(fft2(topoinverse8.^5))))+((((frequencytotal.^(5))/(prod(1:6))).*(fft2(topoinverse8.^6))))+((((frequencytotal.^(6))/(prod(1:7))).*(fft2(topoinverse8.^7))))+((((frequencytotal.^(7))/(prod(1:8))).*(fft2(topoinverse8.^8))))+((((frequencytotal.^(8))/(prod(1:9))).*(fft2(topoinverse8.^9))));

topo9=topo9.*filter;
    						topo9=constant-topo9;
    						topoinverse9=real(ifft2(topo9));
    						diference9=topoinverse9-topoinverse8;
    						diference9=diference9.^2;
      					rms9=sqrt(sum(sum(diference9))/(2*(numrows*numcolumns)));
                     rms=rms9
                     finaltopoinverse=topoinverse9;
                     iter=9;
                     if rms9>=criterio
                        topo10=(((frequencytotal.^(1))/(prod(1:2))).*(fft2(topoinverse9.^2)))+((((frequencytotal.^(2))/(prod(1:3))).*(fft2(topoinverse9.^3))))+((((frequencytotal.^(3))/(prod(1:4))).*(fft2(topoinverse9.^4))))+((((frequencytotal.^(4))/(prod(1:5))).*(fft2(topoinverse9.^5))))+((((frequencytotal.^(5))/(prod(1:6))).*(fft2(topoinverse9.^6))))+((((frequencytotal.^(6))/(prod(1:7))).*(fft2(topoinverse9.^7))))+((((frequencytotal.^(7))/(prod(1:8))).*(fft2(topoinverse9.^8))))+((((frequencytotal.^(8))/(prod(1:9))).*(fft2(topoinverse9.^9))))+((((frequencytotal.^(9))/(prod(1:10))).*(fft2(topoinverse9.^10))));

topo10=topo10.*filter;
    						topo10=constant-topo10;
    						topoinverse10=real(ifft2(topo10));
    						diference10=topoinverse10-topoinverse9;
    						diference10=diference10.^2;
                     rms10=sqrt(sum(sum(diference10))/(2*(numrows*numcolumns)));
                     rms=rms10
                     finaltopoinverse=topoinverse10;
                     iter=10;
                   end
                end
             end
          end
       end
    end
   end
  end

iter     %it displays the iteration number at which the process has stopped
rms      %it displays the rms error at the end of the iterative process 

%the matrix is transposed to their output
trans=(finaltopoinverse)';

%At this point, the forward modelling starts. The final topography obtained before is used in the Parker's formula in order to compute the gravity anomaly due to that interface. The number of terms of the series is the same that the number of iterations reached in the inverse modelling.

sumas=1; %initialize variables. The maximum number of sums is 10, because the maximum number of iterations in the inverse modelling is 10.
sumtotal=[];
sum2=[];
sum3=[];
sum4=[];
sum5=[];
sum6=[];
sum7=[];
sum8=[];
sum9=[];
sum10=[];
constantbouinv=-(2*pi*6.67*contrast)*(exp(-(z0*100000).*(frequencytotal.*(1/100000))));  %the constant term of the series is computed
%The sum of fourier transforms starts here. At each step, a new term of the series is computed and added to the previous one. Then, if the number of sums is lower than the number of iterations, the process continues.

sumtotal=(((frequencytotal.^(0))/(prod(1:1))).*(fft2(finaltopoinverse.^1)));
sums=2;
if sums<=iter
   sum2=(((frequencytotal.^(1))/(prod(1:2))).*(fft2(finaltopoinverse.^2)));
   sums=3;
   sumtotal=sumtotal+sum2;
   if sumas<=iter
      sum3=((((frequencytotal.^(2))/(prod(1:3))).*(fft2(finaltopoinverse.^3))));
      sums=4;
      sumtotal=sumtotal+sum3;
      if sums<=iter
         sum4=((((frequencytotal.^(3))/(prod(1:4))).*(fft2(finaltopoinverse.^4))));
         sums=5;
         sumtotal=sumtotal+sum4;
         if sums<=iter
            sum5=((((frequencytotal.^(4))/(prod(1:5))).*(fft2(finaltopoinverse.^5))));    		
            sums=6;
            sumtotal=sumtotal+sum5;
            if sums<=iter
               sum6=((((frequencytotal.^(5))/(prod(1:6))).*(fft2(finaltopoinverse.^6))));
               sums=7;
               sumtotal=sumtotal+sum6;
               if sums<=iter
                  sum7=((((frequencytotal.^(6))/(prod(1:7))).*(fft2(finaltopoinverse.^7))));
                  sums=8;
                  sumtotal=sumtotal+sum7;
                  if sums<=iter
                     sum8=((((frequencytotal.^(7))/(prod(1:8))).*(fft2(finaltopoinverse.^8))));
                     sums=9;
                     sumtotal=sumtotal+sum8;
                     if sums<=iter
                        sum9=((((frequencytotal.^(8))/(prod(1:9))).*(fft2(finaltopoinverse.^9))));
                        sums=10;
                        sumtotal=sumtotal+sum9;
                        if sums<=iter
                           sum10=((((frequencytotal.^(9))/(prod(1:10))).*(fft2(finaltopoinverse.^10))));
                           sumtotal=sumtotal+sum10;
                        end
                     end
                  end
               end
            end
         end
      end
   end
end


sumtotal=sumtotal.*constantbouinv;  %after the summatory, the final map in frequency domain is computed multipliying the constant term
sumtotal(1,1)=0;
bouinv=real(ifft2(sumtotal));  %the inverse 2-D FFT provides the gravity map in space domain

grdwrite2(Xg,Yg,finaltopoinverse,['Uplate.' num2str(ii) '.grd']);

% reformat grid (output from grdwrite2 is deprecated and does not cooperate
% with grdflexure)
system(['grd2xyz Uplate.' num2str(ii) '.grd | xyz2grd -GUplate.' num2str(ii) '.grd -R' grdfile]);
