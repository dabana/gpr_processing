function [seis,t,x]=plotimage_getseismic
% Returns the seismic matrix in true amplitude and the t and x coordinates
% Create from PI_rescale with some modifications
global SMEAN STDDEV MXS
mainax=findobj(gcf,'type','axes','tag','MAINAXES');
%posax=findobj(gcf,'type','axes','tag','POSITIONAXES');
% try
h=get(gcf,'userdata');
hmsg=h(2);
hdat=h(4);
hi=h(5);
hscale=h(6);
hclip=h(7);
hmaster=h(10);
hlimbox=h(14);
ampflag=get(hmaster,'value');
stuff=get(hmaster,'userdata');
mxs=stuff(1);
smean=stuff(2);
stddev=stuff(3);
%get the data
smat=get(hi,'cdata');
%determine old scaling
dat=get(hscale,'userdata');
oldscaleopt=dat(1);
mxsold=dat(2); mnsold=dat(3); smeanold=dat(4); stddevold=dat(5);
%new opt
newscaleopt=get(hscale,'value');
dat(1)=newscaleopt;
%get clip value
inewclip=get(hclip,'value');
ioldclip=get(hclip,'userdata');
set(hclip,'userdata',inewclip);
clips=get(hdat,'userdata');
clipold=clips(ioldclip);
clipnew=clips(inewclip);

%get number of columns in colormap
clrmap=get(gcf,'colormap');
nkols=size(clrmap,1);

%flag=computer;
%flag='shit';

if(ampflag==1)
    % Figure is Independant
    col=[.8314 .8157 .7843];
    TTS='';
    stddev2=stddev;
    mxs2=mxs;
    smean2=smean;
elseif(ampflag==2)
    % Figure is Master
    col=[1 .5 0];
    TTS='';
    %global SMEAN STDDEV MXS
    stddev2=stddev;
    mxs2=mxs;
    smean2=smean;
    SMEAN=smean2;
    STDDEV=stddev;
    MXS=mxs;
elseif(ampflag==3)
    % while figure is slaved, automatically Max Scaling
    col=[1 1 0];
    TTS='Scaling Automatically Forced Master Figure Scaling';
    % Figure is Slave
    limlns=findobj(gcf,'type','line','tag','LIMITLINE');
    limpts=findobj(gcf,'type','line','tag','LIMITPOINT');
    limdat=get(hlimbox,'userdata');
    if(~isempty(limlns))
        delete(limlns);
        delete(limpts);
        delete(limdat{3});
        set(hlimbox,'userdata',[]);
    end
    %global SMEAN STDDEV MXS
    stddev2=STDDEV;
    mxs2=MXS;
    smean2=SMEAN;
end
set(hmaster,'backgroundcolor',col,'tooltipstring',TTS);
% 		%undo the old scaling
if( oldscaleopt == 1 ) %undo mean scaling
    mxsprime = min([smeanold+clipold*stddevold,mxsold]);
    mns=-mxsprime;
    seis = (smat-1)*(mxsprime-mns)/(nkols-1) + mns;
elseif( oldscaleopt == 2) %undo max scaling
    mns=-mxsold;
    seis = (smat-1)*(mxsold-mns)/(nkols-1) + mns;
end

x=get(hi,'xdata');
t=get(hi,'ydata');

