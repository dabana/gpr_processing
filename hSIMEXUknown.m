function [h,rho,indrho]=hSIMEXUknown(W,Y,errortype,sigU)

%Author: Aurore Delaigle
%Computes bandwidth h and ridge parameter rho using a version of the SIMEX method of
%Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice in errors-in-variables problems. JASA, 103, 280-287 
%
%WARNING: these are not the codes used in the original paper. This is a simplified version of those codes.
%
%Use the function NWDecUknown to compute the regression estimator with this rho and this h

%W: vector of contaminated data W_1,...,W_n
%Y: vector of data Y_1,...,Y_n
%h: bandwidth
%
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%sigU: parameter of Laplace or normal errors used only to define phiU.
%rho: ridge parameter. 


% --------------------------------------------------------
% Preliminary calculations and initialisation of functions
% --------------------------------------------------------


%Default values of phiU(t)=characteristic function of the errors
%If you want to consider another error type, simply replace phiU by the characteristic function of your error type

if strcmp(errortype,'Lap')==1
	phiU=@(t) 1./(1+sigU^2*t.^2);
	varU=2*sigU^2;
elseif strcmp(errortype,'norm')==1
	phiU = @(t) exp(-sigU^2*t.^2/2);
	varU=sigU^2;
end


%phiK: Fourier transform of the kernel K. You can change this if you wish to use another kernel but make sure 
%you change the range of t-values, which should correspond to the support of phiK
phiK = @(t) (1-t.^2).^3;


%Range of t-values (must correspond to the domain of phiK)
dt = .0002;
t = (-1:dt:1)';
t=reshape(t,length(t),1);



n=length(W);
W=reshape(W,1,n);
Y=reshape(Y,1,n);

%number of bins used to compute CV in each SIMEX world
nbin=min(100,n);

%Number of SIMEX samples
BB=20;

%Define a grid where to search for the SIMEX bandwidth. By default we take [h/2,2h], where h=PI bandwidth for density estimation.
%Increase the gird if too small
hPIfX=PI_deconvUknownth4(W,errortype,varU,sigU);
a=hPIfX/2;
b=2*hPIfX;
gridh=a:(b-a)/20:b;


%Define a defaul grid where to search for rho. 
%Recall that rho prevents the denominator of the NW estimator from being too small. In the SIMEX world, the denominator estimates the contaminated density f_W
%This is what motivates the default grid for rho used here.

%Estimator of fW(q_{0.05}) and fW(q_{0.95}) using standard (error-free) KDE and normal reference bandwidth, where q_{alpha} denotes the alpha empirical quantile of the W_i's.
hW=1.06*sqrt(var(W))*n^(-1/5);
ab=quantile(W,[0.05,0.95]);
xout=outerop(ab,W,'-');
fWEF=normpdf(xout,0,hW)*ones(n,1)/n;
gridrho=min(fWEF)*(0.025:0.025:4);


lh=length(gridh);
lrho=length(gridrho);
CVrho=zeros(lh,lrho);
	


%---------------------------------------------------------------------
%Step 1: find the ridge parameter using only the first level of SIMEX
%---------------------------------------------------------------------

%Bin the W data to speed up the computations
[midbin,indbin]=BinData(W,nbin);

for bb=1:BB


	%Generate SIMEX data Wstar
	if strcmp(errortype,'Lap')==1
		Wstar=W+rlap(sigU,1,n);
	elseif strcmp(errortype,'norm')==1
		Wstar=W+normrnd(0,sigU,1,n);
	end

	%For each h in the grid of h-candidates, compute the CV criterion for the data Wstar (this will automatically consider all rho candiates)
	for kh=1:lh
		h=gridh(kh);
		CVrho(kh,:)=CVrho(kh,:)+NWDecridgeL1OCUknown(Wstar,Y,errortype,sigU,h,gridrho,midbin,indbin,nbin);
	end

end


%find which pair of (h,rho) minimizes CV
minCV=min(find(CVrho==min(min(CVrho))));
[indh,indrho]=ind2sub(size(CVrho),minCV)

%Rigdge parameter
rho=gridrho(indrho);

%h from SIMEX level 1
h1=gridh(indh);


%----------------------------------------
%Step 2: Keep rho fixed and find h SIMEX 
%----------------------------------------

CVhstar=0*gridh;
	
for bb=1:BB

	
	%Generate SIMEX data Wstar2
	if strcmp(errortype,'Lap')==1
		Wstar=W+rlap(sigU,1,n);
		Wstar2=Wstar+rlap(sigU,1,n);
	elseif strcmp(errortype,'norm')==1
		Wstar=W+normrnd(0,sigU,1,n);
		Wstar2=Wstar+normrnd(0,sigU,1,n);
	end

	
	%Bin the Wstar data to speed up the computations
	[midbin,indbin]=BinData(Wstar,nbin);

	%Compute CV for each h in the grid, using the ridge parameter rho found above
	for kh=1:lh
		h=gridh(kh);
		CVhstar(kh)=CVhstar(kh)+NWDecridgeL1OCUknown(Wstar2,Y,errortype,sigU,h,rho,midbin,indbin,nbin);
	end

end


indh=find(CVhstar==min(CVhstar))
	

%h from SIMEX level 2
h2=min(gridh(indh));

%Finally deduce h SIMEX
h=h1^2/h2;