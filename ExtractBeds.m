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
imagesc([d(1) d(end)],[t(1) t(end)],trace_trnk,CLIM0)

%(1)Slant stack artificially incoherent data (with randomly reversed
%trace_trnks) and estimate Pn'(x) locally from histograms

mask=ones(size(trace_trnk,1),size(trace_trnk,2))-repmat(2*randint(size(trace_trnk,2),1)',size(trace_trnk,1),1);
trace_trnk_n=mask.*trace_trnk;
[stp_n,tau_n,p_n]=tptran(trace_trnk_n,t,d,-50,50,.1);%1er iteration
%[SlStk_noise,tau_n,p_n]=tptran(noise,t,d,-50,50,.1);%2eme iteration
%We previously make the slant stack for d' before any statitics
[stp,tau,p]=tptran(trace_trnk,t,d,-50,50,0.1);%1er iteration
%[SlStk_signal,tau,p]=tptran(seisclip,t,d,-50,50,0.1);%2eme iteration

figure
imagesc([tau(1) tau(end)],[p(1) p(end)],stp)
%imagesc(stp)
%imagesc(SlStk_noise)
figure
imagesc([tau_n(1) tau_n(end)],[p(1) p(end)],stp_n)
%imagesc(stp_n)
%imagesc(SlStk_signal)


%We now loop over every sample
windowSize=[10 20];
N = size(stp,1) - windowSize(1) + 1;
M = size(stp,2) - windowSize(2) + 1;
E = zeros(N ,M);
E2 = zeros(N ,M);
ZeroMask = zeros(N ,M);
%%
%[XatPdprimmax, PnprimmuMat, PnprimsigmaMat]=deal(zeros(N,M));%2em itération

for i = 1:N
    for j = 1:M
        tstart=tic;
        %cut-out pieces of stp and stp_n
        sample=reshape(stp(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);%première itération
        sample_n=reshape(stp_n(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);%première itération
        %sample=reshape(SlStk_signal(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);%deuxième itération
        %sample_n=reshape(SlStk_noise(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);%deuxième itération
        %toc(tstart)
        
        %estimate Pnprim by Gaussian fit
        [Pnprimmu,Pnprimsigma]=normfit(sample_n);
        %[PnprimmuMat(i,j),PnprimsigmaMat(i,j)]=normfit(sample_n);%deuxième itération, sauve une matrice des paramètres
        %estimate Pnprim by kernel estimation
        %[bandwidth,Pnprimkde,xmeshn,cdf] = kde(sample_n,2^8,-10e6,10e6);
        %clear xmeshn cdf
        %toc(tstart)
        
        %(2, part1) Slant stack [done previously] data and estimate Pd'(x) locally from histograms.
        [bandwidth,Pdprim,xmesh,cdf] = kde(sample,2^8,-10e6,10e6);
        %toc(tstart)
        
        %Creating a mask for the step «zero all samples containing
        %significant percentages of noise»
        noise_threshold=0.25;
        [dump,ind]=max(Pdprim);%1er itération
        clear dump
        %XatPdprimmax(i,j)=xmesh(ind);%2em itération sauve une matrice de x@max(Pdprim)
        %noise_threshold=5e5;
        
        if abs(xmesh(ind)-Pnprimmu)/Pnprimsigma>noise_threshold;
        %end %With this «end» we evaluate deconvolution for all the image

        %(2, part2) Estimate Ps(x) from step (1). [using deconvolution]
        %déconvolution by kernel densitity estimation from Aurore Ladaigle
        hPI = PI_deconvUknownth4(sample,'norm',Pnprimsigma,Pnprimsigma);
        fXdecUK=fdecUknown(xmesh,sample,hPI,'norm',Pnprimsigma,diff(xmesh(1:2)));%bandwith ou hPI?
        %toc(tstart)
        
        %(3 part 1) Evaluate E(s|d') for each sample of the transformed data.
        dprim=stp(i,j);
        [dump,I]=min(abs(dprim-xmesh));
        clear dump
        PdprimATdprim=Pdprim(I);
        Pnprim=normpdf(dprim-xmesh,Pnprimmu,Pnprimsigma);
        E(i,j)=sum(xmesh.*fXdecUK.*Pnprim*diff(xmesh(1:2)))/PdprimATdprim;%première iétration
        %E2(i,j)=sum(xmesh.*fXdecUK.*Pnprim*diff(xmesh(1:2)))/PdprimATdprim;%deuxième iétration
        toc(tstart)
        end %With this end instead of the earlier one, we evaluate the 
        %constly deconvolution only for low noise samples
        

        
        %Visualize distributions
        close all
        figure
        plot(xmesh,Pdprim,'b');
        hold all
        plot(dprim-xmesh,Pnprim,'k');
        %plot(xmesh,Pnprimkde,'k');
        plot(xmesh,fXdecUK,'r');
        plot(E(i,j),PdprimATdprim,'.','markerSize',12);
        plot(dprim,PdprimATdprim,'.','markerSize',12);
        
    end
end
%%
%plot E(s|d')
figure
CLIM=[-8e6 8e6];
imagesc([min(tau_E) max(tau_E)],[min(p_E) max(p_E)],E,CLIM);

%creating the mask
out_n_th=2.5e5;
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

%criterium
Sigma_th0=1.8;
Sigmatau=0.01;
%ZeroMaskAnalytic=0.5+0.5*erf((SignNoiseMeasure2-Sigma_th0)/Sigmatau);%Usign the normalised signal mesure criterium

 %close all
% figure
% plot(reshape(overlap,size(ZeroMask,1)*size(ZeroMask,2),1),reshape(ZeroMaskAnalytic,size(ZeroMask,1)*size(ZeroMask,2),1),'.')
% figure
% plot(sort(reshape(overlap,size(ZeroMask,1)*size(ZeroMask,2),1)))


Emask=E.*ZeroMaskAnalytic;
% figure
% imagesc(Emask)

%smooth an array of E(s|d')/d' values both spatially and temporally before
%multiplying by d' vector
Eond=Emask./stp(1:end-windowSize(1)+1,1:end-windowSize(2)+1);
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
Esmth=smooth2a(Eond,150,25).*stp(1:end-windowSize(1)+1,1:end-windowSize(2)+1);%.*ZeroMaskAnalytic;%it helps to smooth more in p then in tau

%(4 part 1) Invert the extracted signal... 

%inverse slat stack the transformed (exctracted) data
tau_E=tau(1:1:end-windowSize(1)+1);
p_E=p(1:1:end-windowSize(2)+1);
[seis,t1,x1]=itptran(Esmth,tau_E,p_E,d(1),d(end),.2);

%(4 part 2)... and substract from the original data

%Truncate seis both in time and values
Iendoftime=find((t(1)+t1)>t(end),1,'first');
seisclip=seis(1:Iendoftime,:);
%Show seis clipped spatially
% figure
% imagesc([x2(1) x2(end)],[t(1)+t2(1) t(1)+t2(Iendoftime)],seisclip);

for i=1:size(seisclip,1)
    for j=1:size(seisclip,2)
        if abs(seis(i,j))>32767
        seisclip(i,j)=sign(seis(i,j))*32767;
        end
    end
end


%Show seis_clip
figure
clim=[-37000 37000];
imagesc([x2(1) x2(end)],[t(1)+t2(1) t(1)+t2(Iendoftime)],seisclip,clim);

%subtract seis_clip to original data and show the resulting noise
noise=trace_trnk-seisclip;
figure
clim=[-10000 10000];
imagesc([x2(1) x2(end)],[t(1)+t2(1) t(1)+t2(Iendoftime)],noise);





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




