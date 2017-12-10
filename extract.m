function [Signal,ZeroMask,SignNoiseMeasure] = extract(d,n,windowSize,smthMat,out_n_th,in_n_th)
%Optimized function to extrat signal
%OUTPUTS

%Signal - Extracted signal matix
%ZeroMask - Mask used for extraction
%SignNoiseMeasure - Matrix containing the difference between noise and
%signal for every samples. It is used to create ZeroMask.

%INPUTS

%d - Sample matrix
%n - noise estimation matrix
%windowSize - matrix [M N] conting the dimensions o the window over which
%the local statistics are calculated
%out_n_th - outter noise threashold
%in_n_th - inner noise threshold
%smthMat - an array [Vsmth Hsmth] specifying the extend of vertical
%smoothing (Vsmth) and horizontal smoothing (Hsmth).

%-------------------------------------------------------------------------%


%Initialize some variables
N = size(d,1) - windowSize(1) + 1;
M = size(d,2) - windowSize(2) + 1;
E = zeros(N ,M);
SignNoiseMeasure = zeros(N ,M);

for i = 1:N
    for j = 1:M
        tstart=tic;
        %cut-out pieces of d and n
        sample_d=reshape(d(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        sample_n=reshape(n(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        %toc(tstart)
        
        %estimate Pnprim by Gaussian fit
        [Pnmu,Pnsigma]=normfit(sample_n);
        %estimate Pnprim by kernel estimation
        %[bandwidth,Pnprimkde,xmeshn,cdf] = kde(sample_n,2^8,-10e6,10e6);
        %clear xmeshn cdf
        %toc(tstart)
        
        %(2, part1)estimate Pd(x) locally from histograms.
        [bandwidth,Pd,xmesh,cdf] = kde(sample_d,2^8,-10e6,10e6);
        clear bandwidth cdf
        %toc(tstart)
        
        %Update E(i,j) (with costly deconvolution) only for those samples containing few noise as
        %defined by out_n_th. For other samples, E(i,j)=0.
        [dump,ind]=max(Pd);
        clear dump
        SignNoiseMeasure(i,j)=abs(xmesh(ind)-Pnmu);
        if SignNoiseMeasure(i,j)>in_n_th
            
            %(2, part2) Estimate Ps(x) from step (1). [using deconvolution]
            %deconvolution by kernel densitity estimation from Aurore Ladaigle
            hPI = PI_deconvUknownth4(sample_d,'norm',Pnsigma,Pnsigma);
            fXdecUK=fdecUknown(xmesh,sample_d,hPI,'norm',Pnsigma,diff(xmesh(1:2)));%bandwith ou hPI?
            %toc(tstart)
            
            %(3 part 1) Evaluate E(s|d') for each sample of the transformed data.
            [dump,I]=min(abs(d(i,j)-xmesh));
            clear dump
            PdATd=Pd(I);
            Pn=normpdf(d(i,j)-xmesh,Pnmu,Pnsigma);
            E(i,j)=sum(xmesh.*fXdecUK.*Pn*diff(xmesh(1:2)))/PdATd;
            %toc(tstart)
        end
        
    end
end

%Zero all samples containing significant percentages of noise
%creating the zeroing mask
ZeroMask=SignNoiseMeasure>out_n_th;
%smoothing it
h=fspecial('gaussian',[5 5]);
ZeroMaskAnalytic=imfilter(ZeroMask,h);
%applying it
Emask=E.*ZeroMaskAnalytic;

%smooth an array of E(s|d')/d' values both spatially and temporally before
%multiplying by d' vector
Eond=Emask./d(1:(end-windowSize(1)+1),1:(end-windowSize(2)+1));
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
Signal=smooth2a(Eond,smthMat(1),smthMat(2)).*d(1:(end-windowSize(1)+1),1:(end-windowSize(2)+1));


end

