%Script pour extraire les lits horizontaux et inclinés
%d'après Harlan WS et al. 1984

%definition axe distances
dist=0.2*(0:1:2568);
%Selection de la boite d'intéret
itmin=150;
itmax=350;
ixmin=565;
ixmax=865;
trace_trnk=trace(itmin:itmax,ixmin:ixmax);
t=time(itmin:itmax);
d=dist(ixmin:ixmax);

%Show boîte d'intérêt location
figure
CLIM0=[-1000 1000];
imagesc(trace,CLIM0)
hold all
plot([ixmin ixmax],[itmin itmin],'k-','LineWidth',2)
plot([ixmin ixmax],[itmax itmax],'k-','LineWidth',2)
plot([ixmax ixmax],[itmin itmax],'k-','LineWidth',2)
plot([ixmin ixmin],[itmin itmax],'k-','LineWidth',2)
%Show boîte d'intérêt zoom
figure
CLIM0=[-37000 37000];
imagesc([d(1) d(end)],[t(1) t(end)],trace_trnk,CLIM0)


%Define noise and signal
[stp_n,tau_n,p_n]=tptran(noise,t,d,-50,50,.1);
%We previously make the slant stack for d' before any statitics
[stp,tau,p]=tptran(seisclip,t,d,-50,50,0.1);

%(Part 2 to 3) Taken care of by the extract.m function for optimized
%version or functions Es_d follow by ExtractforEsd for non-optimized
%with an adjustable noise threshold

%optimized...
[Signal2,ZeroMask2,SignNoiseMeasure2] = extract(stp,stp_n,[10 20],[100 25],5e5,6.5e5);

%non optimized...

%(4 part 1) Invert the extracted signal... 

%inverse slat stack the transformed (exctracted) data
tau_E=tau(1:end-9);
p_E=p(1:end-19);
[seis2,t1,x1]=itptran(Signal2,tau_E,p_E,d(1),d(end),.2);

%(4 part 2)... and substract from the original data

%Truncate seis both in time and values
Iendoftime=find((t(1)+t1)>t(end),1,'first');
seisclip2=seis2(1:Iendoftime,:);
for i=1:size(seisclip2,1)
    for j=1:size(seisclip2,2)
        if abs(seis2(i,j))>32767
        seisclip2(i,j)=sign(seis2(i,j))*32767;
        end
    end
end
%Show seis_clip
figure
clim=[-10000 10000];
imagesc([x2(1) x2(end)],[t(1)+t2(1) t(1)+t2(Iendoftime)],seisclip2);

%subtract seis_clip to original data and show the resulting noise
noise2=trace_trnk-seisclip2;
figure
clim=[-10000 10000];
imagesc([x2(1) x2(end)],[t(1)+t2(1) t(1)+t2(Iendoftime)],noise2);

%show histogram of the transformed (extracted) data
figure
seis_rh=reshape(seis,size(seis,1)*size(seis,2),1);
hist(seis_rh,100);

%Show histogram of original data
figure
trace_trnk_rh=reshape(trace_trnk,size(trace_trnk,1)*size(trace_trnk,2),1);
hist(trace_trnk_rh,100);

%Show histogram of inversed slant stacked original data
figure
%inverse slant stack the original data
[seis_o,t2,x2]=itptran(stp,tau,p,d(1),d(end),.2);
seis_o_rh=reshape(seis_o,size(seis_o,1)*size(seis_o,2),1);
hist(seis_o_rh,100);

%Show inverse slant stack of original data
figure
clim=[-10000 10000];
imagesc([x2(1) x2(end)],[t(1)+t2(1) t(1)+t2(end)],seis_o);

%Show histogram of exctraction
figure
noise_rh=reshape(noise,size(noise,1)*size(noise,2),1);
hist(noise_rh,100);




