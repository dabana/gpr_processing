function hPI = PI_deconvUknownth4(W,errortype,varU,sigU)

%Author: Aurore Delaigle
%compute 2-stage plug-in bandwidth for kerndel deconvolution estimator as in:
%Delaigle, A. and I. Gijbels (2002). Estimation of integrated squared density derivatives from a contaminated sample, Journal of the Royal Statistical Society, B, 64, 869-886.
%Delaigle, A. and I. Gijbels (2004). Practical bandwidth selection in deconvolution kernel density estimation, Computational Statistics and Data Analysis, 45, 249 - 267

%W: vector of contaminated data
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%varU: variance of errors U_i
%sigU: parameter of Laplace or normal errors used only to define phiU.


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


%Default values of phiU(t)=characteristic function of the errors
%If you want to consider another error type, simply replace phiU by the characteristic function of your error type

if strcmp(errortype,'Lap')==1
	phiU=@(t) 1./(1+sigU^2*t.^2);
elseif strcmp(errortype,'norm')==1
	phiU = @(t) exp(-sigU^2*t.^2/2);
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


%Sample size
n = length(W);

%grid of h values where to search for a solution: you can change the default grid if no solution is found in this grid.
maxh=(max(W)-min(W))/10;

%NR bandwidth of the KDE estimator using the same kernel as we use in the DKDE case
hnaive=((8*sqrt(pi)*RK/3/muK2^2)^0.2)*sqrt(var(W))*n^(-1/5);
hgrid=hnaive/3:(maxh-hnaive/3)/100:maxh;
lh = length(hgrid);
hgrid=reshape(hgrid,1,lh);


%Estimator of the standard deviation of X
stdevx = max(sqrt(var(W) - varU),1/n);



%Quantities that will be needed several times in the computations below
toverh=t*(1./hgrid);
phiK2=(phiK(t)).^2;
phiU2=phiU(toverh).^2;




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
term2=term2./(2*pi*n*hgrid.^(2*rr+1));

ABias2 = (term1 + term2).^2;

%Print the index of the minimiser of Abias2 to see if we are inside the grid (if not, enlarge the grid of bandwidths)
indh3=find(ABias2==min(ABias2),1,'first');
h3 = hgrid(indh3);

%Estimate empirical characteristic function of W
OO=outerop(t/h3,W,'*');
rehatphiW=sum(cos(OO),2)/n;
imhatphiW=sum(sin(OO),2)/n;
clear OO;

%Compute th3
normhatphiW2=rehatphiW.^2+imhatphiW.^2;
th3 = sum(t.^(2*rr) .* normhatphiW2 .* phiK2 ./ phiU2(:,indh3));
th3 = th3*deltat/(2*pi*h3^(2*rr+1));

% -----------------------------------------------------
% Find bandwidth h2 for computing th2, then compute th2
% -----------------------------------------------------


rr=2;
term1= -hgrid.^2*muK2*th3;
term2=repmat(t.^(2*rr).*phiK2,1,lh)./phiU2;
term2=sum(term2,1)*deltat./(2*pi*n*hgrid.^(2*rr+1));

ABias2 = (term1 + term2).^2;

%Print the index of the minimiser of Abias2 to see if we are inside the grid (if not, enlarge the grid of bandwidths)
indh2=find(ABias2==min(ABias2),1,'first');
h2 = hgrid(indh2);

%Estimate empirical characteristic function of W
OO=outerop(t/h2,W,'*');
rehatphiW=sum(cos(OO),2)/n;
imhatphiW=sum(sin(OO),2)/n;
clear OO;


%Compute th2
normhatphiW2=rehatphiW.^2+imhatphiW.^2;
th2 = sum(t.^(2*rr) .* normhatphiW2 .* phiK2 ./ phiU2(:,indh2));
th2=th2*deltat/(2*pi*h2^(2*rr+1));


% ------------------------------------------------------------------------------------------------------
% Finally, compute the bandwidth that minimises the AMISE of the deconvolution kernel density estimator
% ------------------------------------------------------------------------------------------------------

term1=hgrid.^4*muK2^2*th2/4;
term2=repmat(phiK2,1,lh)./phiU2;
term2=sum(term2,1)*deltat./(2*pi*n*hgrid);
AMISE=term1+term2;

indh=find(AMISE==min(AMISE),1,'first');
hPI = hgrid(indh);



