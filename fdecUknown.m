function fXdecUK=fdecUknown(xx,W,h,errortype,sigU,dx)

%Author: Aurore Delaigle
%Compute the deconvolution kernel density estimator

%xx: vector of x-values where to compute the deconvolution kernel density estimator
%W: contaminated sample
%h bandwidth
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%sigU: parameter of Laplace or normal errors used only to define phiU.
%
%
%rescale: to rescale the estimator so that it integrates to 1 after the negative parts have been truncated to zero
%Rescaling requires xx to be a fine grid of equispaced x-values that covers the whole range of x-values where the estimated density is significantly non zero.
%In this case dx=distance between two neighbour points in the xx grid

rescale=0;
if nargin==6
	rescale=1;
end



% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%								WARNINGS:
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% The range of t-values -1 and 1 correspond to the support of phiK. 
% If you change phiK and take a kernel for which phiK is not supported on [-1,1] you have to change -1 and 1 accordingly
%
% muK2 below is the second moment \int x^2 K(x) dx  of the kernel K defined below through phiK. 
% If you change phiK, you MUST change muK2 accordingly.
%
% The phiK here must match the phiK used to compute the bandwidth (PI, CV or other)
%
% The DKDE can also be computed using the Fast Fourier Transform, which is a bit more complex. 
% See Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating integrals and optimizing objective functions: a case study in density deconvolution.   Statistics and Computing,  17,  349 - 355
% However if the grid of t-values is fine enough, the estimator can simply be computed like here without having problems with oscillations.
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------



%Default values of phiU(t)=characteristic function of the errors
%If you want to consider another error type, simply replace phiU by the characteristic function of your error type

if strcmp(errortype,'Lap')==1
	phiU=@(t) 1./(1+sigU^2*t.^2);
elseif strcmp(errortype,'norm')==1
	phiU=@(t) exp(-sigU^2*t.^2/2);
end

%phiK: Fourier transform of the kernel K. You can change this if you wish to use another kernel but make sure 
%you change the range of t-values, which should correspond to the support of phiK
phiK = @(t) (1-t.^2).^3;


%Range of t-values (must correspond to the domain of phiK)
%If you get poor results, check that this grid is fine enough
deltat = .02;%0.0002 is default
t = (-1:deltat:1)';
t=reshape(t,length(t),1);




n=length(W);
OO=outerop(t/h,W,'*');
phiUth=phiU(t/h);



%Estimate empirical characteristic function of W
rehatphiX=sum(cos(OO),2)./phiUth/n;
imhatphiX=sum(sin(OO),2)./phiUth/n;

xt=outerop(t/h,xx,'*');
longx=length(xx);

%Compute the DKDE estimator
fXdecUK=cos(xt).*repmat(rehatphiX,1,longx)+sin(xt).*repmat(imhatphiX,1,longx);
fXdecUK=sum(fXdecUK.*repmat(phiK(t),1,longx),1)/(2*pi)*deltat/h;
fXdecUK(fXdecUK<0)=0*fXdecUK(fXdecUK<0);

if rescale==1
	fXdecUK=fXdecUK/sum(fXdecUK)/dx;
end