
timeMat=repmat(time(66:length(time),:),1,num_trace);
%Creation DES matrices de vitesse
num_trace=size(trace,2); %total number of traces in the profile

%Création de la matrice des vitesses constance
Velocity=ones(1,num_trace)*0.11; %vecteur de vitesse en (m/ns)
VelocityMat=repmat(Velocity,length(time)-65,1); %IMPORTANT: il y a 65 points avec temps négatifs

%DepthMat2: Calcul de la profondeur valide pour toutes profondeurs (si V(z)=cst ou V(z)=Vrms(z))
DepthMat2=real(sqrt((timeMat.*VelocityMat.*0.5).^2-(0.5^2))); %V(z)=cst
%DepthMat2=real(sqrt((timeMat.*Vrms_mat_trnk.*0.5).^2-(0.5^2))); %V(z)=Vrms(z)

%GET TOPOGRAPHY
Coords=0:headers(2,num_trace)/(num_trace-1):headers(2,num_trace);
ftopo=fopen('line3.top');
tempread=fscanf(ftopo,'%f');
Topo(1,:)=tempread(1:2:end,:)';
Topo(2,:)=tempread(2:2:end,:)';
%Interpolation linéaire des données topographiques
Topointerp=interp1(Topo(1,:),Topo(2,:),Coords','linear');
TopIntSmth=smooth(Topointerp,10)';
TopoMat=repmat(TopIntSmth,length(time)-65,1);

%CORRECTION DU PROFIL
%entrer les variables utiles pour construire l'image
resoZ=0.04; %resolution verticale en m/pixel 
resoX=0.2; %résolution latérale en m/pixel (5 scans par metre)
ele0=165; %Élevation de la première ligne de l'image
elef=70; %Élevation de la dernière ligne de l'image
NbPixZ=(ele0-elef)/resoZ;
NbPixX=num_trace*0.2/resoX;
Image=NaN(NbPixZ,NbPixX);
TimeWindow=length(time(66:end));

%Construction de l'image du radiogramme avec Depthmat2
ElevMat=TopoMat-DepthMat2;
  for i=1:NbPixZ
      diff=abs(ElevMat-(ele0-i*resoZ*ones(TimeWindow,NbPixX)));
      logic=diff==repmat(min(diff),TimeWindow,1);
      Image(i,:)=sum(logic.*trace(66:length(time),:),1);
      clear diff logic
  end
  
%Mask outside
diffImage_smth=smooth2a(diff(Image),10,10);
logic2=diffImage_smth~=0;
Image=Image.*[double(logic2);zeros(1,size(Image,2))];

  %display image
  figure
  CLIM = [-37000 37000];
  imagesc(Image,CLIM);
  hold on
  plot(Coords/resoX,(ele0-TopIntSmth)/resoZ);
  
  %Plot contact au roc
  contacts=[];
  contacts=[floor(230.7/resoX);floor(410.2/resoX);floor(455.7/resoX)];
  contacts=[contacts [(ele0-TopoMat(1,contacts(1,1))+25)/resoZ;(ele0-TopoMat(1,contacts(2,1))+21.56)/resoZ;0]];
  plot(contacts(:,1),contacts(:,2),'k.','MarkerSize',35);
  
  %Plot contact Mb/Gs
%   contacts=[contacts;contacts];
%   contacts(4:end,2)=[(ele0-TopoMat(1,contacts(1,1))+13.28)/resoZ;(ele0-TopoMat(1,contacts(2,1))+5.94)/resoZ;(ele0-TopoMat(1,contacts(2,1))+4.94)/resoZ];
%   plot(contacts(:,1),contacts(:,2),'k.','MarkerSize',35);
  
  %Écrire une image avec alpha rendering
amax=max(max(Image));
amin=min(min(Image));
ImageNorm=mat2gray(Image,[amin amax]);
ImageAlpha=mat2gray(abs(Image),[0 1000]);
imwrite(ImageNorm,'LINE03_bg+migr.png','Alpha',ImageAlpha);

%écrire un SEG-Y en depth pour importer dans GOCAD
altwritesegy('LINE03_depth_img.sgy',Image,resoZ,[],[])