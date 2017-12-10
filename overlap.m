overlap2 = zeros(N ,M);
for i = 1:N
    for j = 1:M
        tstart=tic;
        %cut-out pieces of d and n
        %sample_d=reshape(d(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        %sample_n=reshape(n(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
        sample_d=reshape(SlStk_signal(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);%deuxième itération
        sample_n=reshape(SlStk_noise(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);%deuxième itération
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
        overlap2(i,j)=sum(min(Pd,normpdf(xmesh,Pnmu,Pnsigma)')*diff(xmesh(1:2)));
    end
end
