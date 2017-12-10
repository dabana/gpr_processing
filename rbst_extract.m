%function [Signal,overlap] = rbst_extract(d,n,windowSize,smthMat,OL_th,OL_th0,OLtau)
%Optimized function to extrat signal
%OUTPUTS

%Signal - Extracted signal matix
%Overlap - Matrix containing the overlap integrals between Pd and Pn

%INPUTS
%d - Sample matrix
%n - noise estimation matrix
%windowSize - matrix [M N] conting the dimensions o the window over which
%the local statistics are calculated
%OL_th - overlap threshold, overlap area over which Signal is zero and
%costly deconvolution is not computed
%OL_th0 - overlap threshold for the definition of the erf mask
%OLtau - Abruptness of the erf mask
%smthMat - an array [Vsmth Hsmth] specifying the extend of vertical
%smoothing (Vsmth) and horizontal smoothing (Hsmth).

%-------------------------------------------------------------------------%


%Initialize some variables
N = size(d,1) - windowSize(1) + 1;
M = size(d,2) - windowSize(2) + 1;
E = zeros(N ,M);
overlap = zeros(N ,M);

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
        
        %calculate overlap integral
        overlap(i,j)=sum(min(Pd,normpdf(xmesh,Pnmu,Pnsigma)')*diff(xmesh(1:2)));
    end
end

%Update E(i,j) (with costly deconvolution) only for those samples containing few noise as
        %defined by OL_th. For other samples, E(i,j)=0.
        
for i = 1:N
    for j = 1:M
    
        if overlap(i,j)<OL_th
            
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

%applying it
ZeroMask=0.5-0.5*erf((overlap-OL_th0)/OLtau);
Emask=E.*ZeroMask;

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


%end

