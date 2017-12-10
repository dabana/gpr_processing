%creating the signal to noise measure based on abs(x@max(Pdprim)-Pnprimu)
SignNoiseMeasure=zeros(N,M);
for i = 1:N
    for j = 1:M
        sample=reshape(stp(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        sample_n=reshape(stp_n(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        [bandwidth,Pdprim,xmesh,cdf] = kde(sample,2^8,-10e6,10e6);
        
        [Pnprimmu,Pnprimsigma]=normfit(sample_n);
        
        [dump,ind]=max(Pdprim);
        clear dump
        SignNoiseMeasure(i,j)=abs(xmesh(ind)-Pnprimmu);
    end
end

%creating the signal to noise measure based on abs(x@max(Pdprim)-Pnprimu)/Pnprimsigma 
SignNoiseMeasure2_2=zeros(N,M);
for i = 1:N
    for j = 1:M
        %sample=reshape(stp(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        %sample_n=reshape(stp_n(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        sample=reshape(SlStk_signal(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);%deuxième itération
        sample_n=reshape(SlStk_noise(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);%deuxième itération
        [bandwidth,Pdprim,xmesh,cdf] = kde(sample,2^8,-10e6,10e6);
        
        [Pnprimmu,Pnprimsigma]=normfit(sample_n);
        
        [dump,ind]=max(Pdprim);
        clear dump
        SignNoiseMeasure2_2(i,j)=abs(xmesh(ind)-Pnprimmu)/Pnprimsigma;
    end
end

%creating the mask
ZeroMask4=zeros(N,M);
noiseratio_threshold=7.9e5;
for i = 1:N
    for j = 1:M
        if SignNoiseMeasure(i,j)<noiseratio_threshold;
        ZeroMask4(i,j)=0;
        else
        ZeroMask4(i,j)=1;
        end
    end
end

%creating the mask for criterium 2
ZeroMask5=zeros(N,M);

for i = 1:N
    for j = 1:M
        if SignNoiseMeasure2(i,j)<1.75;
        ZeroMask5(i,j)=0;
        else
        ZeroMask5(i,j)=1;
        end
    end
end

%show histogram of the signal to noise measure
figure
SignNoiseMeasure_rh=reshape(SignNoiseMeasure,size(SignNoiseMeasure,1)*size(SignNoiseMeasure,2),1);
hist(SignNoiseMeasure_rh,100);