function [E,SignNoiseMeasure,SignNoiseMeasure2] = Es_d(d,n,windowSize,in_n_th)
%{
Function to calculate E(s|d) the esperance value of the signal s knowing the
samples d and a noise estimate n.

OUTPUTS

E - Extracted E(s|d) matix

SignNoiseMeasure - Matrix containing the absolute value of the difference between noise and
signal for every samples (abs(xmesh(ind)-Pnmu)). It is used to create the
ZeroMask in ExctractforEsd.m.

SignNoiseMeasure2 - Matrix containing, for every samples, the absolute value of the difference between noise and
signal normalized by the standard deviation of the noise (abs(xmesh(ind)-Pnmu)/Pnsigma). It CAN be used to create the
ZeroMask in ExctractforEsd.m but is less effective than SignNoiseMeasure.
Not that to recover he standard deviation of the noise, you just do
SignNoiseMeasure./SignNoiseMeasure2.

INPUTS

d - Sample matrix

n - noise estimation matrix

windowSize - matrix [M N] conting the dimensions o the window over which
the local statistics are calculated

in_n_th - inner noise threashold defined as a multiple of the standard
deviation of noise. If the maximum of the sample PDF is located within
mean +- in_n_th*sigma, the deconvolution is not calculated.

smthMat - an array [Vsmth Hsmth] specifying the extend of vertical
smoothing (Vsmth) and horizontal smoothing (Hsmth).
%}
%-------------------------------------------------------------------------%


%Initialize some variables
N = size(d,1) - windowSize(1) + 1;
M = size(d,2) - windowSize(2) + 1;
E = zeros(N ,M);
SignNoiseMeasure = zeros(N ,M);
SignNoiseMeasure2 = zeros(N ,M);
ranged=[min(min(d)) max(max(d))];

for i = 1:N
    for j = 1:M
        %tstart= tic;
        %cut-out pieces of d and n
        sample_d=reshape(d(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        sample_n=reshape(n(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        
        %(1 part 2)... Estimate Pn(x) locally from histograms (or by Gaussian fit in our
        %case)
        [Pnmu,Pnsigma]=normfit(sample_n);
        
        %(2, part1) estimate Pd(x) localy from histograms (or kernel
        %density estimation in our case)
        [~,Pd,xmesh,~] = kde(sample_d,2^8,-10e6,10e6);
                
        %Calculation deconvolution only for samples containing
        %significantly less noise based on SignNoiseMeasure2
        [~,ind]=max(Pd);
        SignNoiseMeasure(i,j)=abs(xmesh(ind)-Pnmu);
        SignNoiseMeasure2(i,j)=abs(xmesh(ind)-Pnmu)/Pnsigma;
        
        if SignNoiseMeasure2(i,j)>in_n_th
            %(2, part2) Estimate Ps(x) from step (1). [using deconvolution]
            %deconvolution by kernel densitity estimation from Aurore Ladaigle
            hPI = PI_deconvUknownth4(sample_d,'norm',Pnsigma,Pnsigma);
            fXdecUK=fdecUknown(xmesh,sample_d,hPI,'norm',Pnsigma,diff(xmesh(1:2)));%bandwith ou hPI?

            %(3 part 1) Evaluate E(s|d') for each sample of the transformed data.
            [~,I]=min(abs(d(i,j)-xmesh));
            PdATd=Pd(I);
            Pn=normpdf(d(i,j)-xmesh,Pnmu,Pnsigma);
            E(i,j)=sum(xmesh.*fXdecUK.*Pn*diff(xmesh(1:2)))/PdATd;

        end
        
        
        
        %toc(tstart)
        
    end
    disp([num2str(i) '/' num2str(N)])
    disp([num2str(j) '/' num2str(M)])
   
end
end