
%Plot location of hyperbolas over profil
%close all
plot(Vrms(:,2)/.2,52.472+Vrms(:,3)/.8,'.');
text(Vrms(:,2)/.2,52.472+Vrms(:,3)/.8,num2str(Vrms(:,4)),'VerticalAlignment','bottom','HorizontalAlignment','right')

%Open the GRD file with krigging from velocity in a matrix
%Don't forget saving as Surfer 6 text grid
[Vrms_matrix xmin xmax ymin ymax]=grd_read_v2('L0+00fortes_GPRv.grd');
logic=Vrms_matrix~=max(max(Vrms_matrix));
Vrms_matrix=[Vrms_matrix.*logic;zeros(1125-size(Vrms_matrix,1),2569)];
figure
imagesc(Vrms_matrix);
hold all
plot(Vrms1(:,2)/.2,52.472+Vrms1(:,3)/.8,'.');
Vrms1=Vrms(Vrms(:,1)==1,:);
text(Vrms1(:,2)/.2,52.472+Vrms1(:,3)/.8,num2str(Vrms1(:,4)),'VerticalAlignment','bottom','HorizontalAlignment','right')

%Creation DES matrices de vitesse
num_trace=size(trace,2); %total number of traces in the profile

%Création de la matrice des vitesses constance
Velocity=ones(1,num_trace)*0.11; %vecteur de vitesse en (m/ns)
VelocityMat=repmat(Velocity,length(time)-65,1); %IMPORTANT: il y a 65 points avec temps négatifs

%CALCUL DES MATRICE D'ÉLEVATION
%Calcul d'une matrice d'élevation en approche simple
Vrms_mat_trnk=Vrms_matrix(66:length(time),:);
timeMat=repmat(time(66:length(time),:),1,num_trace);
%DepthMat=0.5*timeMat.*Vrms_mat_trnk; %IMPORTANT timeMat = two-way time!
DepthMat=0.5*timeMat.*VelocityMat; %IMPORTANT timeMat = two-way time!

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
  
  figure
  CLIM = [-1000 1000];
  imagesc(Image,CLIM);
  %Plot contact au roc
  contacts=[];
  contacts=[floor(230.7/resoX);floor(410.2/resoX);floor(455.7/resoX)];
  contacts=[contacts [(ele0-TopoMat(1,contacts(1,1))+25)/resoZ;(ele0-TopoMat(1,contacts(2,1))+21.56)/resoZ;0]];
  
  %Plot contact Mb/Gs
  contacts=[contacts;contacts];
  contacts(4:end,2)=[(ele0-TopoMat(1,contacts(1,1))+13.28)/resoZ;(ele0-TopoMat(1,contacts(2,1))+5.94)/resoZ;(ele0-TopoMat(1,contacts(2,1))+4.94)/resoZ];
  hold all
  plot(contacts(:,1),contacts(:,2),'k.','MarkerSize',35);
  
  %Construction de l'image du model de vitesses avec Depthmat
ElevMat=TopoMat-DepthMat;
  for i=1:NbPixZ
      diff=abs(ElevMat-(ele0-i*resoZ*ones(TimeWindow,NbPixX)));
      logic=diff==repmat(min(diff),TimeWindow,1);
      ImageVrms(i,:)=sum(logic.*Vrms_matrix(66:length(time),:),1);
      clear diff logic
  end
   figure
   imagesc(ImageVrms);
     hold all
  plot(contacts(:,1),contacts(:,2),'k.','MarkerSize',35);

  
