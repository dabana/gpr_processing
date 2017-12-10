function CV=NWDecridgeL1OCUknown(W,Y,errortype,sigU,h,rhogrid,midbin,indbin,nbin)

%Author: Aurore Delaigle
	%Compute a version of weighted CV used in SIMEX, for binned data
	%W plays the role of the contaminated data
	%midbin is the vector of the centers of the binned data that play the role of the non contaminated data

	%Default values of phiU(t)=characteristic function of the errors
	%If you want to consider another error type, simply replace phiU by the characteristic function of your error type

	if strcmp(errortype,'Lap')==1
		phiU=@(t) 1./(1+sigU^2*t.^2);
	elseif strcmp(errortype,'norm')==1
		phiU = @(t) exp(-sigU^2*t.^2/2);
	end


	%phiK: Fourier transform of the kernel K. You can change this if you wish to use another kernel but make sure 
	%you change the range of t-values, which should correspond to the support of phiK		
	%The phiK used in this code must be the same as in the other codes (hSIMEX + NW codes)
	phiK = @(t) (1-t.^2).^3;
	dt = .001;
	t = (-1:dt:1)';
	th=t/h;
	longt=length(t);
	

	n=length(W);

	
	%Compute the empirical characteristic function of W (times n) at t/h: \hat\phi_W(t/h)
	OO=(ones(n,1)*th(:).');
	OO=W(:)*ones(1,longt).*OO;
	csO=cos(OO);
	snO=sin(OO);
	clear OO;
	
	rehatphiW=sum(csO,1);
	imhatphiW=sum(snO,1);

	%Compute \sum_j Y_j e^{itW_j/h}
	Y=reshape(Y,1,n);
	renum=Y*csO;
	imnum=Y*snO;

	%Compute \hat m(M_i) where M_i is the middle of the bin in which X_i (the non contaminated data) lies
	xt=th(:)*ones(1,nbin);     
	xt=xt.*(ones(longt,1)*midbin(:).');
	cxt=cos(xt);
	sxt=sin(xt);
	clear xt;
	cxt=cxt(:,indbin);
	sxt=sxt(:,indbin);
	

	phiUth=phiU(th);
	matphiKU=reshape(phiK(t)./phiUth,1,longt);
	Den=(rehatphiW.*matphiKU)*cxt+(imhatphiW.*matphiKU)*sxt;
	Num=(renum.*matphiKU)*cxt+(imnum.*matphiKU)*sxt;


	%Compute from there the leave-one-out version \hat m_{-i}(M_i)
	
	csO=csO';
	snO=snO';
	
	Den=Den-matphiKU*(csO.*cxt)-matphiKU*(snO.*sxt);
	for i=1:n
		csO(:,i)=csO(:,i)*Y(i);
		snO(:,i)=snO(:,i)*Y(i);
	end
	
	Num=Num-matphiKU*(csO.*cxt)-matphiKU*(snO.*sxt);


	%Finally compute weighted CV, where the ith term of the sum is weighted by f_W(W_i)
	
	rhogrid=rhogrid*(2*pi*h*n)/dt;
	hW=1.06*sqrt(var(W))*n^(-1/5);
	xout=outerop(W,W,'-');
	fWEF=normpdf(xout,0,hW)*ones(n,1)/n;
	
	CV=0*rhogrid;
	for krho=1:length(rhogrid)
		rho=rhogrid(krho);
		dd=Den;
		dd(dd<rho)=rho;

		mhatstar=Num./dd;
		mhatstar=reshape(mhatstar,1,n);		
		CV(krho)=sum(fWEF'.*(Y-mhatstar).^2);
	end

