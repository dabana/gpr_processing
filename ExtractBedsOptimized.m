%Script pour extraire les lits horizontaux et inclinés
%d'après Harlan WS et al. 1984
close all
%definition axe distances
dist=0.2*(0:1:2568);
%Selection de la boite d'intéret
%{
itmin=1;
itmax=length(time);
ixmin=870;
ixmax=930;
%}
itmin=90;
itmax=150;
ixmin=2000;
ixmax=2150;
trace_trnk=trace(itmin:itmax,ixmin:ixmax);
t=time(itmin:itmax);
x=dist(ixmin:ixmax);

%Show boîte d'intérêt location
figure
CLIM0=[-40000 40000];
imagesc(trace,CLIM0)
hold all
plot([ixmin ixmax],[itmin itmin],'k-','LineWidth',2)
plot([ixmin ixmax],[itmax itmax],'k-','LineWidth',2)
plot([ixmax ixmax],[itmin itmax],'k-','LineWidth',2)
plot([ixmin ixmin],[itmin itmax],'k-','LineWidth',2)
plot(V_emp(:,1)/0.2,(V_emp(:,2)+52.472)/0.8,'.k','MarkerSize',20)
text(V_emp(:,1)/0.2,(V_emp(:,2)+52.472)/0.8,num2str(V_emp(:,3)), 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right','color',[1 1 1]);
%Show boîte d'intérêt zoom
figure
CLIM0=[-37000 37000];
imagesc([x(1) x(end)],[t(1) t(end)],trace_trnk,CLIM0)
hold all
plot(V_emp(:,1),V_emp(:,2),'.k','MarkerSize',20)
text(V_emp(:,1),V_emp(:,2),num2str(V_emp(:,3)), 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right','color',[1 1 1]);
%%

IT=2;
%(1 part 1) Slant stack artificially incoherent data... (with randomly reversed
%trace_trnks)


if IT==1
mask=ones(size(trace_trnk,1),size(trace_trnk,2))-repmat(2*randint(size(trace_trnk,2),1)',size(trace_trnk,1),1);
trace_trnk_n=mask.*trace_trnk;
[stp_n,tau_n,p_n]=tptran(trace_trnk_n,t,x,-10,10,0.015);
[stp,tau,p]=tptran(trace_trnk,t,x,-10,10,0.015);
else
[stp,tau,p]=tptran(trace_trnk,t,x,-10,10,0.015);
[stp_n,tau_n,p_n]=tptran(noise,t,x,-10,10,0.015);
end

%We previously make the slant stack for d' before any statistics
%[stp,tau,p]=tptran(trace_trnk,t,x,-10,10,0.025);


%Vizualise slant stack data and noise
figure
CLIM0=[-10e5 10e5];
imagesc([p(1) p(end)],[tau(1) tau(end)],stp,CLIM0)
figure
CLIM0=[-10e5 10e5];
imagesc([p_n(1) p_n(end)],[tau_n(1) tau_n(end)],stp_n,CLIM0)

%invert for fun
%[seis,t1,x1]=itptran(stp,tau,p,x(1),x(end),.2);

%%
%(Steps 1 part 2, 2 and 3) Taken care of by Es_d.m followed by ExtractforEsd.m
windowSize=[30 5];

%[E,SignNoiseMeasure,SignNoiseMeasure2] = Es_d(stp,stp_n,windowSize,.5);
%Here is a piece of Es_d
%Initialize some variables
d=stp;
n=stp_n;
N = size(d,1) - windowSize(1) + 1;
M = size(d,2) - windowSize(2) + 1;
[Pnmu, Pnsigma, SignNoiseMeasure, SignNoiseMeasure2, E] = deal(zeros(N ,M));
Pd=cell(N ,M);
kde_min=1.1*min(min(stp));
kde_max=1.1*max(max(stp));

for i = 1:N
    for j = 1:M
        %tstart= tic;
        %cut-out pieces of d and n
        sample_d=reshape(d(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        sample_n=reshape(n(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        
        %(1 part 2)... Estimate Pn(x) locally from histograms (or by Gaussian fit in our
        %case)
        [Pnmu(i,j),Pnsigma(i,j)]=normfit(sample_n);
        
        %(2, part1) estimate Pd(x) localy from histograms (or kernel
        %density estimation in our case)
        [~,Pd{i,j},xmesh,~] = kde(sample_d,2^8,kde_min,kde_max);
                
        %Calculation deconvolution only for samples containing
        %significantly less noise based on SignNoiseMeasure2
        [~,ind]=max(Pd{i,j});
        SignNoiseMeasure(i,j)=abs(xmesh(ind)-Pnmu(i,j));
        SignNoiseMeasure2(i,j)=abs(xmesh(ind)-Pnmu(i,j))/Pnsigma(i,j);
    end
        disp([num2str(i) '/' num2str(N)])
    disp([num2str(j) '/' num2str(M)])
end
rshp_SNM=reshape(SignNoiseMeasure,M*N,1);

%vizualize signal noise measure
figure
imagesc([p(1) p(end)],[tau(1) tau(end)],SignNoiseMeasure)
%Histogram of signal noise measure
[Nbin,Xbin]=hist(rshp_SNM,500);
out_n_th=Xbin(find(sum(repmat(Nbin',1,500).*tril(ones(500,500)))<(0.05*N*M),1,'first'));

%DECONVOLUTION
%vizualize outter mask
figure
imagesc([p(1) p(end)],[tau(1) tau(end)],SignNoiseMeasure>out_n_th)
[Ikeep,Jkeep]=find(SignNoiseMeasure>out_n_th);
for k=1:length(Ikeep)
    
            tstart=tic;
            sample_d=reshape(d(Ikeep(k):Ikeep(k)+windowSize(1)-1,Jkeep(k):Jkeep(k)+windowSize(2)-1),windowSize(1)*windowSize(2),1);
            
           
            %(2, part2) Estimate Ps(x) from step (1). [using deconvolution]
            %deconvolution by kernel densitity estimation from Aurore Ladaigle
            hPI = PI_deconvUknownth4(sample_d,'norm',Pnsigma(Ikeep(k),Jkeep(k)),Pnsigma(Ikeep(k),Jkeep(k)));
            fXdecUK=fdecUknown(xmesh,sample_d,hPI,'norm',Pnsigma(Ikeep(k),Jkeep(k)),diff(xmesh(1:2)));%bandwith ou hPI?

            %(3 part 1) Evaluate E(s|d') for each sample of the transformed data.
            [~,I]=min(abs(d(Ikeep(k),Jkeep(k))-xmesh));
            PdATd=Pd{Ikeep(k),Jkeep(k)}(I);
            Pn=normpdf(d(Ikeep(k),Jkeep(k))-xmesh,Pnmu(Ikeep(k),Jkeep(k)),Pnsigma(Ikeep(k),Jkeep(k)));
            E(Ikeep(k),Jkeep(k))=sum(xmesh.*fXdecUK.*Pn*diff(xmesh(1:2)))/PdATd; 
            
            toc(tstart)
end



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
%ExtractforEsd(E,stp(1:end-windowSize(1)+1,1:end-windowSize(2)+1),SignNoiseMeasure,7.9e5,[100
%25]); %Cette fonction là ne marche pas bien...
%Signal=E;

%creating the mask
in_n_th=2.75*out_n_th;
ZeroMask=zeros(size(E,1),size(E,2));
for i = 1:size(E,1)
    for j = 1:size(E,2)
        if SignNoiseMeasure(i,j)<in_n_th;
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
Signal=smooth2a(Eond,100,25).*stp(1:end-windowSize(1)+1,1:end-windowSize(2)+1);%.*ZeroMaskAnalytic;%it helps to smooth more in p then in tau

%(4 part 1) Invert the extracted signal... 

%inverse slat stack the transformed (exctracted) data
tau_E=tau(1:end-windowSize(1)+1);
p_E=p(1:end-windowSize(2)+1);
[seis,t1,x1]=itptran(Signal,tau_E,p_E,x(1),x(end),.2);

%(4 part 2)... and substract from the original data

%Truncate seis both in time and values
Iendoftime=find((t(1)+t1)>t(end),1,'first');
seisclip=seis(1:Iendoftime-1,:);
treal=t(1)+t1(1:Iendoftime-1);

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

% close all
% factor=[0.5 0.75 1 1.25 1.5 2];
% for h=1:length(factor)
%     noise=trace_trnk-factor(h)*seisclip;
%     CorNoiseTrace{h}=xcorr2(noise,trace_trnk);
%     CorNoiseSeis{h}=xcorr2(noise,seisclip);
% end

% figure
% for h=1:length(factor)
% subplot(2,3,h)
% imagesc(CorNoiseTrace{h})
% title(sprintf(['factor= ' num2str(factor(h))]))
% end
% 
% figure
% for h=1:length(factor)
% subplot(2,3,h)
% imagesc(CorNoiseSeis{h})
% title(sprintf(['factor= ' num2str(factor(h))]))
% end
% 
% close all
% factor=[0.5 0.75 1 1.25 1.5 2];
% for h=1:length(factor)
%     noise=trace_trnk-factor(h)*seisclip;
%     this=corrcoef(noise(:),trace_trnk(:));
%     CorCoefNoiseTrace(h)=this(1,2);
%     this=corrcoef(noise(:),seisclip(:));
%     CorCoefNoiseSeis(h)=this(1,2);
% end
% figure
% plot(factor,CorCoefNoiseSeis);
% hold all
% plot(factor,-CorCoefNoiseTrace);
% figure
% plot(factor,abs(CorCoefNoiseTrace)./CorCoefNoiseSeis);
% 


close all
%Show seis_clip
figure
clim=[-10000 10000];
imagesc([x1(1) x1(end)],[t(1)+t1(1) t(1)+t1(Iendoftime)],seisclip);

%subtract seis_clip to original data and show the resulting noise
figure
clim=[-10000 10000];

C=0.25;
if IT==2
noise2=trace_trnk-C*seisclip;
imagesc([x1(1) x1(end)],[t(1)+t1(1) t(1)+t1(Iendoftime)],noise2);
else
noise=trace_trnk-C*seisclip;
imagesc([x1(1) x1(end)],[t(1)+t1(1) t(1)+t1(Iendoftime)],noise);
end

%Show boîte d'intérêt zoom
figure
CLIM0=[-37000 37000];
imagesc([x(1) x(end)],[t(1) t(end)],trace_trnk,CLIM0)
