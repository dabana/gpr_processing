%Script pour extraire les lits horizontaux et inclinés
%d'après Harlan WS et al. 1984

%definition axe distances
dist=0.2*(0:1:2568);
%Selection de la boite d'intéret
%{
itmin=1;
itmax=length(time);
ixmin=870;
ixmax=930;
%}
itmin=150;
itmax=400;
ixmin=675;
ixmax=775;
trace_trnk=trace(itmin:itmax,ixmin:ixmax);
t=time(itmin:itmax);
x=dist(ixmin:ixmax);

%Show boîte d'intérêt location
figure
CLIM0=[-37000 37000];
imagesc(trace,CLIM0)
hold all
plot([ixmin ixmax],[itmin itmin],'k-','LineWidth',2)
plot([ixmin ixmax],[itmax itmax],'k-','LineWidth',2)
plot([ixmax ixmax],[itmin itmax],'k-','LineWidth',2)
plot([ixmin ixmin],[itmin itmax],'k-','LineWidth',2)
%Show boîte d'intérêt zoom
figure
CLIM0=[-37000 37000];
imagesc([x(1) x(end)],[t(1) t(end)],trace_trnk,CLIM0)
%%
%(1 part 1) Slant stack noise

[stp_n,tau_n,p_n]=tptran(noise,t,x,-70,70,0.1);

%We previously make the slant stack for d' before any statistics
[stp,tau,p]=tptran(seisclip,t,x,-70,70,0.1);

%%
%(Steps 1 part 2, 2 and 3) Taken care of by Es_d.m followed by ExtractforEsd.m
WindowSize=[10 20];
[E,SignNoiseMeasure,SignNoiseMeasure2] = Es_d(stp,stp_n,WindowSize,.5);
%get rid of high values?
% N = size(E,1);
% M = size(E,2);
% E2=zeros(N ,M);
% ranged=[min(min(stp)) max(max(stp))];
% for i = 1:N
%     for j = 1:M
% if E(i,j)<ranged(1) || E(i,j)>ranged(2)
%          E2(i,j)=0;
% else
%     E2(i,j)=E(i,j);
% end
%     end
% end
%%
%Signal =
%ExtractforEsd(E,stp(1:end-WindowSize(1)+1,1:end-WindowSize(2)+1),SignNoiseMeasure,7.9e5,[100
%25]); %Cette fonction là ne marche pas bien...
%Signal=E;

%creating the mask
out_n_th=6.9e5;
ZeroMask=zeros(size(E,1),size(E,2));
for i = 1:size(E,1)
    for j = 1:size(E,2)
        if SignNoiseMeasure(i,j)<out_n_th;
        ZeroMask(i,j)=0;
        else
        ZeroMask(i,j)=1;
        end
    end
end

%zero all samples containing significant percentages of noise
h=fspecial('gaussian',[5 5]);
ZeroMaskAnalytic=imfilter(ZeroMask,h);
Emask=E.*ZeroMaskAnalytic;

%smooth an array of E(s|d')/d' values both spatially and temporally before
%multiplying by d' vector
Eond=Emask./stp(1:end-WindowSize(1)+1,1:end-WindowSize(2)+1);
%applying it

[Inan,Jnan] = find(isnan(Eond));
[Iinf,Jinf] = find(isinf(Eond));
[Ibig,Jbig] = find(abs(Eond)>100);
Iout=[Inan;Iinf;Ibig];
Jout=[Jnan;Jinf;Jbig];
for i=1:size(Iout,1)
    for j=1:size(Jout,1)
        Eond(Iout(i),Jout(j))=0;
    end
end
Signal=smooth2a(Eond,100,25).*stp(1:end-WindowSize(1)+1,1:end-WindowSize(2)+1);%.*ZeroMaskAnalytic;%it helps to smooth more in p then in tau

%%
%(4 part 1) Invert the extracted signal... 

%inverse slat stack the transformed (exctracted) data
tau_E=tau(1:end-WindowSize(1)+1);
p_E=p(1:end-WindowSize(2)+1);
[seis,t1,x1]=itptran(Signal,tau_E,p_E,x(1),x(end),.2);

%(4 part 2)... and substract from the original data

%Truncate seis both in time and values
Iendoftime=find((t(1)+t1)>t(end),1,'first');
seisclip=seis(1:Iendoftime,:);

if max(max(abs(seis)))>32767
    %Truncate seis in value
for i=1:size(seisclip,1)
    for j=1:size(seisclip,2)
        if abs(seis(i,j))>32767
        seisclip(i,j)=sign(seis(i,j))*32767;
        end
    end
end
else%Crank up seis
factor=32767/max(max(abs(seis)));
seisclip=seisclip.*factor;
end
%
close all
%Show seis_clip
figure
clim=[-10000 10000];
imagesc([x1(1) x1(end)],[t(1)+t1(1) t(1)+t1(Iendoftime)],seisclip);
figure
clim=[-10000 10000];
imagesc([x1(1) x1(end)],[t(1)+t1(1) t(1)+t1(Iendoftime)],noise);


%subtract seis_clip to original data and show the resulting noise
noise2=noise+seisclip;
figure
clim=[-10000 10000];
imagesc([x1(1) x1(end)],[t(1)+t1(1) t(1)+t1(Iendoftime)],noise2);

%Show boîte d'intérêt zoom
figure
CLIM0=[-37000 37000];
imagesc([x(1) x(end)],[t(1) t(end)],trace_trnk,CLIM0)
