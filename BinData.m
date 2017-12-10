function [midbin,indice]=BinData(W,nbin)

%Author: Aurore Delaigle
%This program bins the data W into nbins
n=length(W);
W=reshape(W,n,1);

%Compute extremities of the bins
ab=quantile(W,[0.05,0.95]);


%Bin widths
delta=(ab(2)-ab(1))/nbin;

%Bin centers
midbin=ab(1)+delta/2+delta*(0:(nbin-1));

%Bin left and right extremities
Abin=midbin-delta/2;
Bbin=midbin+delta/2;


%Find in which bin each observation lies
Wmat=repmat(W,1,nbin);
Amat=repmat(Abin,n,1);
Bmat=repmat(Bbin,n,1);
indice=repmat(1:nbin,n,1);
indice=sum(indice.*((Wmat>Amat)&(Wmat<=Bmat)),2);

%Put those beyond the extremities at the extremities
indice(W<=Abin(1))=1;
indice(W>Bbin(nbin))=nbin;
