function hPI = PI_deconvUknownth4het(W,errortype,sigUj)

%Author: Aurore Delaigle
%compute 2-stage plug-in bandwidth for heteroscedastic kerndel deconvolution estimator as in:
%Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic error. Bernoulli, 14, 562-579

%W: vector of contaminated data
%errortype: 'Lap' for the case where the error densities are all Laplace densities and 'norm' for the case where the error densities are all normal. 
%For other error distributions, simply redefine phiU below 
%sigUj: vector of length n which contains the parameters of each of the n Laplace or normal errors. Used only to define phiU.


% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%								WARNINGS:
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% The range of t-values -1 and 1 correspond to the support of phiK.  
% If you change phiK and take a kernel for which phiK is not supported on [-1,1] you have to change -1 and 1 accordingly
%
% muK2 below is the second moment \int x^2 K(x) dx  of the kernel K defined below through phiK
% If you change phiK, you MUST change muK2 accordingly.
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------







% --------------------------------------------------------
% Preliminary calculations and initialisation of functions
% --------------------------------------------------------

%Sample size
n = length(W);


%Default values of phiU(t)=characteristic function of the errors
%If you want to consider another error type, simply replace phiU by the characteristic function of your error type

if strcmp(errortype,'Lap')==1
	phiUk=@(t,k) 1./(1+sigUj(k)^2*t.^2);
	varX=mean(W.^2)-(mean(W)).^2-2*sum(sigUj.^2)/n;
elseif strcmp(errortype,'norm')==1
	phiUk=@(t,k) exp(-sigUj(k)^2*t.^2/2);
	varX=mean(W.^2)-(mean(W)).^2-sum(sigUj.^2)/n;
end

%phiK: Fourier transform of the kernel K. You can change this if you wish to use another kernel but make sure 
%you change the range of t-values, which should correspond to the support of phiK
phiK = @(t) (1-t.^2).^3;


%second moment \int x^2 K(x) dx  of the kernel K
muK2 = 6;
RK=1024/3003/pi;

%Range of t-values (must correspond to the domain of phiK)
deltat = .0002;
t = (-1:deltat:1);
t=reshape(t,length(t),1);



%grid of h values where to search for a solution: you can change the default grid if no solution is found in this grid.
maxh=(max(W)-min(W))/10;

%NR bandwidth of the KDE estimator using the same kernel as we use in the DKDE case
hnaive=((8*sqrt(pi)*RK/3/muK2^2)^0.2)*sqrt(var(W))*n^(-1/5);
hgrid=hnaive/3:(maxh-hnaive/3)/100:maxh;
lh = length(hgrid);
hgrid=reshape(hgrid,1,lh);


%Estimator of the standard deviation of X
stdevx = max(sqrt(varX),1/n);



%Quantities that will be needed several times in the computations below
toverh=t*(1./hgrid);
phiK2=(phiK(t)).^2;
matphiU=zeros(length(t),length(hgrid));
phiU2=matphiU;
for k=1:n
	matphiU=phiUk(toverh,k);
	phiU2=phiU2+matphiU.^2;
end




% --------------------------------------------
% Estimate theta4 by normal reference method     
% --------------------------------------------

th4 = stdevx^(-9)*105/(32*sqrt(pi)); 

% ------------------------------------------------------
% Find bandwidth h3 for computing th3, then compute th3
% ------------------------------------------------------

rr=3;
term1= -hgrid.^2*muK2*th4;
term2=repmat(t.^(2*rr).*phiK2,1,lh)./phiU2;
term2=sum(term2,1)*deltat;
term2=term2./(2*pi*hgrid.^(2*rr+1));

ABias2 = (term1 + term2).^2;

%Print the index of the minimiser of Abias2 to see if we are inside the grid (if not, enlarge the grid of bandwidths)
indh3=find(ABias2==min(ABias2),1,'first')
h3 = hgrid(indh3);

%Estimate empirical characteristic function of W
OO=outerop(t/h3,W,'*');

% Compute phiU(-t/h) -- since phiU is symmetric, this is the same as phiU(t/h)
matphiU=OO;
for k=1:n
	matphiU(:,k)=phiUk(t/h3,k);
end

phiUth=sum(matphiU.^2,2);

%Estimate empirical characteristic function of W.
rehatphiX=sum(cos(OO).*matphiU,2)./phiUth;
imhatphiX=sum(sin(OO).*matphiU,2)./phiUth;



clear OO;

%Compute th3
normhatphiX2=rehatphiX.^2+imhatphiX.^2;
th3 = sum(t.^(2*rr) .* normhatphiX2 .* phiK2);
th3 = th3*deltat/(2*pi*h3^(2*rr+1));

% -----------------------------------------------------
% Find bandwidth h2 for computing th2, then compute th2
% -----------------------------------------------------


rr=2;
term1= -hgrid.^2*muK2*th3;
term2=repmat(t.^(2*rr).*phiK2,1,lh)./phiU2;
term2=sum(term2,1)*deltat./(2*pi*hgrid.^(2*rr+1));

ABias2 = (term1 + term2).^2;

%Print the index of the minimiser of Abias2 to see if we are inside the grid (if not, enlarge the grid of bandwidths)
indh2=find(ABias2==min(ABias2),1,'first')
h2 = hgrid(indh2);

%Estimate empirical characteristic function of W
OO=outerop(t/h2,W,'*');
% Compute phiU(-t/h) -- since phiU is symmetric, this is the same as phiU(t/h)
matphiU=OO;
for k=1:n
	matphiU(:,k)=phiUk(t/h2,k);
end

phiUth=sum(matphiU.^2,2);

%Estimate empirical characteristic function of W.
rehatphiX=sum(cos(OO).*matphiU,2)./phiUth;
imhatphiX=sum(sin(OO).*matphiU,2)./phiUth;


clear OO;


%Compute th2
normhatphiX2=rehatphiX.^2+imhatphiX.^2;
th2 = sum(t.^(2*rr) .* normhatphiX2 .* phiK2);
th2=th2*deltat/(2*pi*h2^(2*rr+1));


% ------------------------------------------------------------------------------------------------------
% Finally, compute the bandwidth that minimises the AMISE of the deconvolution kernel density estimator
% ------------------------------------------------------------------------------------------------------

term1=hgrid.^4*muK2^2*th2/4;
term2=repmat(phiK2,1,lh)./phiU2;
term2=sum(term2,1)*deltat./(2*pi*hgrid);
AMISE=term1+term2;

indh=find(AMISE==min(AMISE),1,'first')
hPI = hgrid(indh);



