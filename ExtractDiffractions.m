close all
%if no beds were removed
% noise2=trace_trnk(1:end-windowSize(1)+1,1:end-windowSize(2)+1);
% treal=t(1:end-windowSize(1)+1);
% x1=x(1:end-windowSize(2)+1);
%noise2=noise;

%clip problematic areas
% noise2=noise2(1:end-15,:);
% treal=treal(1:end-15);

%create artificially incoherent data out of noise to create noise_n
mask=ones(size(noise2,1),size(noise2,2))-repmat(2*randint(size(noise2,2),1)',size(noise2,1),1);
noise_n=mask.*noise2;


%Initialize a bunch of variables
v=[0.05e9 0.06e9 0.07e9 0.08e9 0.09e9 0.1e9 0.11e9 0.12e9 0.13e9 0.14e9 0.15e9];%velocities in m/s
[Signal,noise_mig,tmig,xmig,signal_imig,timig,ximig,E_k,SNM_k]=deal(cell(1,length(v)));
windowSize=[15 15];
    
%loop over velocities
for k=1:length(v)
    
    %v(k)=0.15e9;
    %(1) Migrate the data, without continuous bed reflections, over a
    %physical range of velocities
    [noise_mig{k},tmig{k},xmig{k}]=fkmig(noise2,treal.*1e-9,x1,v(k));
    [noise_n_mig,tmig_n,xmig_n]=fkmig(noise_n,treal.*1e-9,x1,v(k));


%Initialize some variables
d=noise_mig{k};
%Check le range de d
n=noise_n_mig;


N = size(d,1) - windowSize(1) + 1;
M = size(d,2) - windowSize(2) + 1;
[Pnmu, Pnsigma, SignNoiseMeasure, E] = deal(zeros(N ,M));
Pd=cell(N ,M);

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
        [~,Pd{i,j},xmesh,~] = kde(sample_d,2^8,-1e5,1e5);
                
        %Calculation deconvolution only for samples containing
        %significantly less noise based on SignNoiseMeasure2
        [~,ind]=max(Pd{i,j});
        SignNoiseMeasure(i,j)=abs(xmesh(ind)-Pnmu(i,j));
        SignNoiseMeasure2(i,j)=abs(xmesh(ind)-Pnmu(i,j))/Pnsigma(i,j);
    end
    disp([num2str(i) '/' num2str(N)])
    disp([num2str(j) '/' num2str(M)])
end
%SNM_k{k}=SignNoiseMeasure;
SNM_k{k}=SignNoiseMeasure2;

rshp_SNM=reshape(SignNoiseMeasure,M*N,1);

%vizualize signal noise measure
% figure
% imagesc([p(1) p(end)],[tau(1) tau(end)],SignNoiseMeasure)
%Histogram of signal noise measure
[Nbin,Xbin]=hist(rshp_SNM,500);
out_n_th=Xbin(find(sum(repmat(Nbin',1,500).*tril(ones(500,500)))<(0.05*N*M),1,'first'));

%DECONVOLUTION
%vizualize outter mask
% figure
% imagesc([p(1) p(end)],[tau(1) tau(end)],SignNoiseMeasure>out_n_th)
[Ikeep,Jkeep]=find(SignNoiseMeasure>out_n_th);

for kk=1:length(Ikeep)
    
            tstart=tic;
            sample_d=reshape(d(Ikeep(kk):Ikeep(kk)+windowSize(1)-1,Jkeep(kk):Jkeep(kk)+windowSize(2)-1),windowSize(1)*windowSize(2),1);
            
           
            %(2, part2) Estimate Ps(x) from step (1). [using deconvolution]
            %deconvolution by kernel densitity estimation from Aurore Ladaigle
            hPI = PI_deconvUknownth4(sample_d,'norm',Pnsigma(Ikeep(kk),Jkeep(kk)),Pnsigma(Ikeep(kk),Jkeep(kk)));
            fXdecUK=fdecUknown(xmesh,sample_d,hPI,'norm',Pnsigma(Ikeep(kk),Jkeep(kk)),diff(xmesh(1:2)));%bandwith ou hPI?

            %(3 part 1) Evaluate E(s|d') for each sample of the transformed data.
            [~,I]=min(abs(d(Ikeep(kk),Jkeep(kk))-xmesh));
            PdATd=Pd{Ikeep(kk),Jkeep(kk)}(I);
            Pn=normpdf(d(Ikeep(kk),Jkeep(kk))-xmesh,Pnmu(Ikeep(kk),Jkeep(kk)),Pnsigma(Ikeep(kk),Jkeep(kk)));
            E(Ikeep(kk),Jkeep(kk))=sum(xmesh.*fXdecUK.*Pn*diff(xmesh(1:2)))/PdATd; 
            
            toc(tstart)
end
E_k{k}=E;
clear Pd Pnmu Pnsigma SignNoiseMeasure
clc
end

%%
%(2) For each migrated section, smoothly zero samples containing
%significant noise 
in_n_th=1*out_n_th;
smthMat=[5 5];
for k=1:length(v)
 
d=noise_mig{k};
%Signal{k} = ExtractforEsd(E_k{k},SNM_k{k},in_n_th,noise_mig{k}(1:(end-windowSize(1)+1),1:(end-windowSize(2)+1)),smthMat);
%Cette fonction là ne marche pas bien...
%Signal=E;

%creating the mask

ZeroMask=zeros(size(E_k{k},1),size(E_k{k},2));
for i = 1:size(E_k{k},1)
    for j = 1:size(E_k{k},2)
        if SNM_k{k}(i,j)<in_n_th;
        ZeroMask(i,j)=0;
        else
        ZeroMask(i,j)=1;
        end
    end
end

%zero all samples containing significant percentages of noise
h=fspecial('gaussian',[5 5]);
ZeroMaskAnalytic=imfilter(ZeroMask,h);
Emask=E_k{k}.*ZeroMaskAnalytic;

%smooth an array of E(s|d')/d' values both spatially and temporally before
%multiplying by d' vector
Eond=Emask./d(1:end-windowSize(1)+1,1:end-windowSize(2)+1);
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
Signal{k}=smooth2a(Eond,smthMat(1),smthMat(2)).*d(1:end-windowSize(1)+1,1:end-windowSize(2)+1);%.*ZeroMaskAnalytic;%it helps to smooth more in p then in tau

end

%%
%(3) Diffract (invert the migration of) each section at the extraction
%velocity

for k=1:length(v)
    
    [signal_imig{k},timig{k},ximig{k}]=ifkmig(Signal{k},tmig{k}(1:(end-windowSize(1)+1)),xmig{k}(1:(end-windowSize(2)+1)),v(k));
%     figure
%     imagesc(signal_imig{k});
end


figure
for k=1:length(v)
subplot(3,4,k)
imagesc(signal_imig{k})
title(sprintf(['v=' num2str(v(k)*1e-9) 'm/ns)']))
end
subplot(3,4,k+1)
imagesc(noise2)


%%
%(4)Find the least-square superposition of these diffracted sections best
%resembling the data without bed reflections 


limit_up=1;
limit_down=size(noise2,1)-10;
matrix=zeros(length(v),length(v));
%Select only the bottom part of the image with parabolaes
noise_clip=noise2(limit_up:(limit_down-windowSize(1)+1),1:(end-windowSize(2)+1));%just the bottom part
MN=size(noise_clip,1)*size(noise_clip,2);


%IMPLEMENTATION OF EQUATION 4
dhat=reshape(noise_clip,MN,1);
source=zeros(1,length(v));

for i=1:length(v)
    e_v=reshape(signal_imig{i}(limit_up:limit_down-windowSize(1)+1,:),MN,1);
    
    source(i)=dot(e_v,dhat);
    for j=1:length(v)
    e_w=reshape(signal_imig{j}(limit_up:limit_down-windowSize(1)+1,:),MN,1);
   
    matrix(i,j)=dot(e_v,e_w);
    end
end

a_w=inv(matrix)*source';

s_hat=zeros(size(noise_clip,1),size(noise_clip,2));
for i=1:length(v)
   s_hat=s_hat+a_w(i).*signal_imig{i}(limit_up:limit_down-windowSize(1)+1,:);
end

%Display the best superposition
figure
imagesc(s_hat,[-50000 50000])
figure
plot(a_w)

%%
%IMPLEMENTATION OF EQUATION E-1
n=4; %order of the polynomial estimate +1
%(3) Diffract (invert the migration of) each section at the extraction
%velocity but be multiplying by t^n before
Gain=zeros(n*length(v),n*length(v));

for k=1:length(v)
    for l=1:n
    ni=(k-1)*n+l;
    %ni=(l-1)*length(v)+k;
    Signal{ni}=E{k};
    Gain(k,l)=((150*tmig{k}(2))^(l-1))^-1;
    [signal_imig{ni},timig{ni},ximig{ni}]=ifkmig(repmat(Gain(k,l).*tmig{k}(1:(end-windowSize(1)+1)).^(l-1),1,282).*Signal{ni},...
        tmig{k}(1:(end-windowSize(1)+1)),xmig{k}(1:(end-windowSize(2)+1)),v(k));
%     figure
%     imagesc(signal_imig{k});
    end
end
limit_up=1;
limit_down=size(noise,1)-0;

matrix=zeros(n*length(v),n*length(v));
%Select only the bottom part of the image with parabolaes
noise_clip=noise(limit_up:(limit_down-windowSize(1)+1),1:(end-windowSize(2)+1));%just the bottom part
MN=size(noise_clip,1)*size(noise_clip,2);
dhat=repmat(reshape(noise_clip,MN,1),1,1);
source=zeros(1,n*length(v));

for ni=1:n*length(v)
    e_v=reshape(signal_imig{ni}(limit_up:limit_down-windowSize(1)+1,:),MN,1);
    source(ni)=dot(e_v,dhat);
    for nj=1:n*length(v)
        e_w=reshape(signal_imig{nj}(limit_up:limit_down-windowSize(1)+1,:),MN,1);
        matrix(ni,nj)=dot(e_v,e_w);
    end
end

close all
imagesc(matrix)
cond(matrix)

a_w=pcg(matrix,source',1e-4,200);
s_hat=zeros(size(noise_clip,1),size(noise_clip,2));
for i=1:n*length(v)
   s_hat=s_hat+a_w(i).*signal_imig{i}(limit_up:limit_down-windowSize(1)+1,:);
end

%Display the best superposition
figure
imagesc(s_hat,[-50000 50000])
%Display the polynomiale weighting as a function of time for a certain
%velocity
v_ind=10;
v_t=Gain(v_ind,1)*a_w(v_ind*n+1)+Gain(v_ind,2)*a_w(v_ind*n+2)*timig{v_ind}+...
        Gain(v_ind,3)*a_w(v_ind*n+3)*timig{v_ind}.^2+Gain(v_ind,4)*a_w(v_ind*n+4)*timig{v_ind}.^3;
figure
plot(timig{1},v_t)
%%
%(5) Migrate this superposition over the previous range of velocities
%(6) Determine the best migration velocity by evaluating the focusing
%measure
Fpx=zeros(length(v),1);
for k=1:length(v)
   
   [s_hat_mig{k},tshatmig_n,xshatmig_n]=fkmig(s_hat,0.8e-9,0.2,v(k));
   %Fpx=zeros(size(s_hat_mig{k},1),size(s_hat_mig{k},2));
   %for i = 1:size(s_hat_mig{k},1)
    %for j = 1:size(s_hat_mig{k},2)
        %estimate P(x) locally from kernel density estimation.

            [~,Px,xmesh,~] = kde(s_hat_mig{k},2^8,-10e4,10e4);
            %[~,Px,xmesh,~] = kde(s_hat_mig{k},2^8);
            Px=abs(Px');
            I=find(Px==0);
            while isempty(I)==0
                I=find(Px==0);
                for l=1:length(I)
                    if I(l)==1
                        Px(I(l))=Px(I(l)+1);
                        
                    elseif I(l)==length(Px)
                        Px(I(l))=Px(I(l)-1);
                        
                    else
                        Px(I(l))=0.5*(Px(I(l)-1)+Px(I(l)+1));
                        
                    end
                    
                end
                I=find(Px==0);
            end
            
            clear bandwidth cdf
            Fpx(k)=sum(Px.*log10(Px)*diff(xmesh(1:2)))+0.5*log10(sum((xmesh.^2).*Px*diff(xmesh(1:2))))+log10(sqrt(2*pi))+0.5;
    %end
   %end
end

%Display the best velocity
[dump,indice]=max(Fpx);
v_best=v(indice)
figure
plot(v,Fpx,'k');
title(sprintf(['Mesure statistique de la focalisation par migration pour trouver V_{optimale}.\nV_{optimale} maximise la mesure de focalisation.\n(ici V_{optimale}=' num2str(v_best*1e-9) 'm/ns)']))
xlabel(sprintf('Vitesse de propagation testée (m/s)'))
ylabel(sprintf('Mesure statistique de la focalisation\npar migration (sans unitées)'))


figure
for k=1:length(v)
subplot(3,4,k)
imagesc(s_hat_mig{k})
title(sprintf(['F[p(x)]= ' num2str(Fpx(k)) 'and v_{best}=' num2str(v(k)*1e-9) 'm/ns)']))
end
subplot(3,4,k+1)
imagesc(noise2)


% %subtract seis_clip to original data and show the resulting noise
%     close all
% noise=trace_trnk-seisclip;
figure
clim=[-10000 10000];
imagesc(noise2);    
% 
% 
%     [seismig,tmig1,xmig1]=fkmig(noise,0.8e-9,.2,.1e9);
%     figure
%     clim=[-10000 10000];
%     imagesc([xmig1(1) xmig1(end)],[tmig1(1) tmig1(end)]*1e9,seismig);
%     
%     [seisimig,timig,ximig]=ifkmig(seismig,tmig1,xmig1,0.1e9);
%     figure
%     clim=[-10000 10000];
%     imagesc([ximig(1) ximig(end)],[timig(1) timig(end)]*1e9,seisimig);
    
    % entropy(k)=sum(sum(seismig.^4,1)./sum(seismig.^2,1),2);