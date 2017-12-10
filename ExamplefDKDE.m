%Put here the path where you have put the Matlab codes, for example:

addpath('C:\Matlab')

%Author: Aurore Delaigle
%This code illustrates how to use the functions for computing the deconvolution kernel density estimator and its bandwidths

%Noise to signal ratio=varU/varX
NSR=0.2

%Sample size
n=500

%Generate data from a normal mixture
X=normrnd(5,.4,1,n);
X2=normrnd(2,1,n);

pmix=0.75;
tmp=unifrnd(0,1,1,n);
X(tmp<pmix)=X2(tmp<pmix);

%true density of X
truedens=@(xx) (1-pmix)*normpdf(xx,5,.4)+pmix*normpdf(xx,2,1);


%Specify error distribution (normal or Laplace in this case) and generate data from this error distribution
errortype='Lap';

if strcmp(errortype,'norm')==1
	%normal case	
	sigU=sqrt(NSR*var(X));
	U=normrnd(0,sigU,1,n);
	varU=sigU^2;


elseif strcmp(errortype,'Lap')==1
	%Laplace case
	sigU=sqrt(NSR*var(X)/2);
	varU=2*sigU^2;
	U=rlap(sigU,1,n);
	
end


%Contaminated data
W=X+U;
varW=var(W);

%Grid where to estimate the true density
xx=-2:0.1:8;
dx=xx(2)-xx(1);

%Plot the true density
plot(xx,truedens(xx),'r')


%PI bandwidth of Delaigle and Gijbels
hPI=PI_deconvUknownth4(W,errortype,varU,sigU);

%DKDE estimator without rescaling (density does not integrate exactly to 1)
y=fdecUknown(xx,W,hPI,errortype,sigU);

%With rescaling: here xx must be equispaced and must cover the range where the estimated density is significantly non zero
y2=fdecUknown(xx,W,hPI,errortype,sigU,dx);

hold
plot(xx,y,'k')
hold

hold
plot(xx,y2,'g')
hold

%CV bandwidth of Stefanski and Carroll
hCV=CVdeconv(W,errortype,sigU)
y3=fdecUknown(xx,W,hCV,errortype,sigU,dx);

hold
plot(xx,y3,'m')
hold


%Compare with the naive KDE estimator that ignores the error (using normal reference bandwidth and standard normal kernel)
h=1.06*sqrt(var(W))*n^(-1/5);
xout=outerop(xx,W,'-');

fnaive=normpdf(xout,0,h)*ones(n,1)/n;
hold
plot(xx,fnaive,'c')
hold

legend('true f', 'fdec, hPI', 'fdec rescaled, hPI', 'fdec rescaled, hCV', 'naive estimator, hNR')

