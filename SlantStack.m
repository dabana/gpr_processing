%close all
%definition axe distances
dist=0.2*(0:1:2568);
%Selection de la boite d'intéret
itmin=150;
itmax=480;
ixmin=20;
ixmax=250;

figure
CLIM0=[-10000 10000];
imagesc(trace,CLIM0)
hold all
plot([ixmin ixmax],[itmin itmin],'k-','LineWidth',2)
plot([ixmin ixmax],[itmax itmax],'k-','LineWidth',2)
plot([ixmax ixmax],[itmin itmax],'k-','LineWidth',2)
plot([ixmin ixmin],[itmin itmax],'k-','LineWidth',2)


trace_trnk=trace(itmin:itmax,ixmin:ixmax);
x=dist(ixmin:ixmax);
t=time(itmin:itmax);
figure
CLIM0=[-10000 10000];
imagesc([x(1) x(end)],[t(1) t(end)],trace_trnk,CLIM0)

%foward slant stack
[stp,tau,p]=tptran(trace_trnk,t,x,-10,10,.1);
figure
imagesc([min(tau) max(tau)],[min(p) max(p)],stp);

%filter out main reflectors by power
threshold=.5e6;
rshp_stp=reshape(stp, prod(size(stp)),1);
rshp_stp=sort(rshp_stp);
stp_filt=stp.*(0.5+0.5.*erf((abs(stp)-threshold)./100));

%statistic plots
figure
imagesc([min(tau) max(tau)],[min(p) max(p)],stp_filt);
figure
plot(rshp_stp)
hold all
plot([0 length(rshp_stp)],[threshold threshold])
figure
hist(rshp_stp,100)
hold all
plot(threshold*[1 1],[0 1e4])

%inverse slant stack of filtered data
[seis,t1,x1]=itptran(stp_filt,tau,p,x(1),x(end),.2);
seis=seis.*repmat((t(1)+t1)<t(end),1,size(seis,2));
figure
CLIM=[-100000 100000];
imagesc([x1(1) x1(end)],[t(1)+t1(1) t(1)+t1(end)],seis,CLIM);

%inverse slant stack of unfiltered data
[seis,t2,x2]=itptran(stp,tau,p,x(1),x(end),.2);
seis=seis.*repmat((t(1)+t1)<t(end),1,size(seis,2));
figure
CLIM=[-100000 100000];
imagesc([x2(1) x2(end)],[t(1)+t2(1) t(1)+t2(end)],seis,CLIM);