%Put here the path where you have put the Matlab codes, for example:

addpath('C:\Matlab')

%Author: Aurore Delaigle
%This code illustrates how to use the functions for computing the deconvolution kernel density estimator and its bandwidths
%in the heteroscedastic case

%Sample size
n=200

mu1=-3;
mu2=2;
sig1=1;
sig2=1;

%Generate data from a normal mixture
X=normrnd(mu1,sig1,1,n);
X2=normrnd(mu2,sig2,1,n);

pmix=0.75;
tmp=unifrnd(0,1,1,n);
X(tmp<pmix)=X2(tmp<pmix);

%true density of X
truedens=@(xx) (1-pmix)*normpdf(xx,mu1,sig1)+pmix*normpdf(xx,mu2,sig2);




%Specify error distributions (normal or Laplace in this case) and generate data from these error distributions
errortype='norm';

sigU=0.6;
%Parameters of the Laplace distribution (one parameter oer observation, so sigmaj is of length n)
sigmaj=sigU*sqrt(1.0+(1:n)/n)*sqrt(0.5);


if strcmp(errortype,'norm')==1
	%normal case	
	U=zeros(1,n);
	for i=1:n
		U(i)=normrnd(0,sigmaj(i),1,1);
	end

elseif strcmp(errortype,'Lap')==1
	%Laplace case
	U=zeros(1,n);
	for i=1:n
	 	U(i)=rlap(sigmaj(i),1,1);
	end
	
end


%Contaminated data
W=X+U;
varW=var(W);

%Grid where to estimate the true density
xx=-8:0.1:7;
dx=xx(2)-xx(1);

%Plot the true density
plot(xx,truedens(xx),'r')


%PI bandwidth of Delaigle and Gijbels
hPI=PI_deconvUknownth4het(W,errortype,sigmaj);

%DKDE estimator without rescaling (density does not integrate exactly to 1)
y=fdecUknownhet(xx,W,hPI,errortype,sigmaj);

%With rescaling: here xx must be equispaced and must cover the range where the estimated density is significantly non zero
y2=fdecUknownhet(xx,W,hPI,errortype,sigmaj,dx);

hold
plot(xx,y,'k')
hold

hold
plot(xx,y2,'g')
hold

%Compare with the naive KDE estimator that ignores the error (using normal reference bandwidth and standard normal kernel)
h=1.06*sqrt(var(W))*n^(-1/5);
xout=outerop(xx,W,'-');

fnaive=normpdf(xout,0,h)*ones(n,1)/n;
hold
plot(xx,fnaive,'c')
hold

legend('true f', 'fdec, hPI', 'fdec rescaled, hPI', 'naive estimator, hNR')

