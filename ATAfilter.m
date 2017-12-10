%ATA filter
ATA=mean(abs(trace)')';
width=50; %number of trace averaged (multiple de 2)
num_trace=size(trace,2); %total number of traces in the profile

%ROLLING AVERAGE IN THE X DIRECTION
% for i=1+width/2:num_trace-width/2
% ATAmatrix(:,i-width/2)=mean(abs(trace(:,i-width/2:i+width/2))')';
% end

%Smoothing ATAmatrix
% smoothw=25;
% for i=1:1:size(ATAmatrix,1)/smoothw-1
%     ATAsmooth(i,:)=mean(ATAmatrix(floor(smoothw*(i-1/2)):floor(smoothw*(i+1/2)),:));
% end

%Calculating energy average trace localy
logATAmatrix=log10(smooth2a(trace.^2,50,250)); %Smoothing first

%calculate attenuation coefficients
% for i=1:size(logATAmatrix,1)/Attw-1
%     tempmat=logATAmatrix(floor(Attw*(i-1/2)):floor(Attw*(i+1/2)),:);
%     for j=1:size(logATAmatrix,2)
%         locfit=fit(time(floor(Attw*(i-1/2)):floor(Attw*(i+1/2))),tempmat(:,j),'poly1');
%         coeff=coeffvalues(locfit);
%         AttlocA(i,j)=coeff(1);
%         AttlocB(i,j)=coeff(2);
%         clear coeff
%         clear locfit
%     end
% end

%DEPTH CALCULATION

%Création de la matrice des vitesses
% Velocity=ones(1,num_trace)*0.110; %vecteur de vitesse en X en (m/ns)
% VelocityMat=repmat(Velocity,length(time)-65,1); %IMPORTANT: il y a 65 points avec temps négatifs

%Création de la matrice des vitesses RMS des hyperboles
% tt=(0.15/0.214)^(-1/0.162);
% VelocityZ=[ones(76,1)*0.150 ; sqrt((((0.15^2)*tt/2)+0.05412.*((time(77:end)./2).^0.676-(tt/2)^0.676))./(time(77:end)./2))];%vecteur de vitesse en Z en (m/ns)
% VelocityMat=repmat(VelocityZ(66:end),1,num_trace);%matrice des vitesse RMS calibrées

%Création de la matrice des vitesses interpolée
load('L000_interface');

%express distance as index
Gs(:,1)=Gs(:,1)./0.2;
GxT(:,1)=GxT(:,1)./0.2;
roc(:,1)=roc(:,1)./0.2;
top=[GxT;Gs];
%Interpolate interfaces on Xindex
Xindex=1:1:num_trace;
top_int=interp1(top(:,1),top(:,2),Xindex,'linear','extrap');
bot_int=interp1(roc(:,1),roc(:,2),Xindex,'linear','extrap');

Velocity=ones(1,num_trace)*0.110;
VelocityMat=repmat(Velocity,length(time),1);
t0=1:1:length(time);

% for k=1:240/0.2
%     tau=(0.126-0.119)/(top_int(k)-bot_int(k));
%     VelocityMat(:,k)=(0.126+(t0-top_int(k))*tau)';    
% end
% 
% for k=(240/0.2)+1:(420/0.2)
%     
%     v1=0.126-(k-(240/0.2))*(0.119-0.126)/((240/0.2)-(420/0.2));
%     v2=0.119-(k-(240/0.2))*(0.115-0.119)/((240/0.2)-(420/0.2));
%         
%     tau=(v1-v2)/(top_int(k)-bot_int(k));
%     VelocityMat(:,k)=(v1+(t0-top_int(k))*tau)';    
% end
% 
% for k=(420/0.2)+1:num_trace
%     tau=(0.119-0.115)/(top_int(k)-bot_int(k));
%     VelocityMat(:,k)=(0.119+(t0-top_int(k))*tau)';    
% end
% 
% VelocityMat=VelocityMat(66:end,:);
% imagesc(VelocityMat);

%Version simplifié
for k=1:240/0.2
    timeindex=round(top_int(k)/0.8); 
    VelocityMat(1:timeindex,k)=0.126;
    VelocityMat(timeindex+1:end,k)=0.119;
end

for k=(240/0.2)+1:(420/0.2)
    
    v1=0.126-(k-(240/0.2))*(0.119-0.126)/((240/0.2)-(420/0.2));
    v2=0.119-(k-(240/0.2))*(0.115-0.119)/((240/0.2)-(420/0.2));
    
    timeindex=round(top_int(k)/0.8);    
    VelocityMat(1:timeindex,k)=v1;
    VelocityMat(timeindex+1:end,k)=v2;  
end

for k=(420/0.2)+1:num_trace
    timeindex=round(top_int(k)/0.8); 
    VelocityMat(1:timeindex,k)=0.119;
    VelocityMat(timeindex+1:end,k)=0.115;    
end

VelocityMat=VelocityMat(66:end,:);
imagesc(VelocityMat);

%Depthmat: Calcul de la profondeur valide où profondeur > séparation des antennes (1
%m)
SampleMat=0.8*ones(length(time)-65,num_trace);%0.8 nanosecondes d'échantillonage
DepthStepmat=0.5.*cat(1,VelocityMat.*SampleMat,zeros(num_trace-length(time)+65,num_trace)); %matrices des pas de profondeurs
%Depthmat=fliplr(rot90(rot90(DepthStepmat)*fliplr(triu(ones(num_trace,length(time)-65)))));

%DepthMat2: Calcul de la profondeur valide pour toutes profondeurs (si V(z)=cst)
%T2wMat=repmat(time(12:end-65),1,num_trace);
%DepthMat2=sqrt((T2wMat.*VelocityMat(12:end,:).*0.5).^2-(0.5^2));
T2wMat=repmat(time(66:end),1,num_trace);
DepthMat2=real(sqrt((T2wMat.*VelocityMat.*0.5).^2-(0.5^2)));

%Calcul de la profondeur hybride de DepthMat2 et Depthmat
%imagesc(100*abs(DepthMat2-Depthmat)./Depthmat); %Erreur entre DepthMat2 et Depthmat
%[row,col]=find(VelocityMat(12:end-1,:)-VelocityMat(13:end,:)); %trouver row la profondeur de la première couche uniforme
clear col

%Depthmat3: Avec un genre hybride (qui ne marche pas car accumule l'error de DepthMat2)
%depthStepMat1stLayer=DepthMat2(2:row,:)-DepthMat2(1:row-1,:);
%DepthStepmat2=[depthStepMat1stLayer;DepthStepmat(row:end,:)];
% Depthmat3=fliplr(rot90(rot90(DepthStepmat2)*fliplr(triu(ones(num_trace,num_trace)))));
% Depthmat3=Depthmat3(12:length(time)-65,:); %11 valeurs au début et 65 à la fin à rejetter, donc 76 au total

%Depthmat4: En collant l'un et l'autre
%Depthmat4=[DepthMat2(1:row,:);Depthmat(row+1:end,:)];

%GET TOPOGRAPHY

ftopo=fopen('line3.top');
tempread=fscanf(ftopo,'%f');
Topo(1,:)=tempread(1:2:end,:)';
Topo(2,:)=tempread(2:2:end,:)';
Coords=Topo(1,1):(headers(2,num_trace)-Topo(1,1))/(num_trace-1):headers(2,num_trace);

%Interpolation linéaire des données topographiques
Topointerp=interp1(Topo(1,:),Topo(2,:),Coords','linear');
TopIntSmth=smooth(Topointerp,50)';
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
      Image(i,:)=sum(logic.*trace(66:end,:),1);
      clear diff logic
  end
  
%   %Construction de l'image du profil smoothé
%   Image2=NaN(NbPixZ,NbPixX);
%    for i=1:NbPixZ
%       di ff=abs(ElevMat-(ele0-i*resoZ*ones(TimeWindow,NbPixX)));
%       logic=diff==repmat(min(diff),TimeWindow,1);
%       Image2(i,:)=sum(logic.*logATAmatrix(12:length(time)-65,:),1);
%       clear diff logic
%   end
 
%Normalisation en image 8-bit
amax=max(max(Image));
amin=min(min(Image));

ImageNorm=mat2gray(Image,[amin amax]);
ImageAlpha=mat2gray(abs(Image),[0 1000]);
figure
imagesc(Image,[-32767,32767]);

%imwrite(ImageNorm,'LINE01.png','Alpha',ImageAlpha);

%exportation vers Surfer!
%grd_write(trace,0,num_trace*resoX,min(time),max(time),'line03_timedomain');%sans correction de profondeur
% grd_write(flipud(Image),-255,273,elef,ele0,'line03.grd');
% grd_write(flipud(Image),-240,288,elef,ele0,'line03_2.grd');

%grd_write(flipud(Image),-245,283,elef,ele0,'line03_v110.grd');
%grd_write(flipud(Image),-245,283,elef,ele0,'line03_Vrms.grd');
grd_write(flipud(Image),-245,283,elef,ele0,'line03_Vmodsimple.grd');

%grd_write(flipud(Image),-100,-100+resoX*num_trace,elef,ele0,'line01_v110.grd');
%grd_write(flipud(Image),-100,-100+resoX*num_trace,elef,ele0,'line01_Vrms.grd');
%grd_write(flipud(Image),210,210+resoX*num_trace,elef,ele0,'line02_v110.grd');
%grd_write(flipud(Image),210,210+resoX*num_trace,elef,ele0,'line02_Vrms.grd');

%grd_write(flipud(Image),0,800,elef,ele0,'line00_Vrms.grd');
%grd_write(flipud(Image),0,800,elef,ele0,'line00_v110.grd');

%grd_write(flipud(Image),0,num_trace*0.2,elef,ele0,'gridline00_Vrms.grd');
%grd_write(flipud(Image),0,num_trace*0.2,elef,ele0,'gridline00_v110.grd');

%grd_write(flipud(Image2),-640,-112,elef,ele0,'line03_smoothed');

  
%%
%Construction de l'image avec ElevMat
%  for i=2:num_trace
%      for j=2:TimeWindow
%      Image(floor((ElevMat(1,1)-ElevMat(j,i))/0.02),i)=trace(j+76,i);
%      end
%  end

%  
%  [X,Y]=meshgrid(1:NbPixZ,1);
% for k=2:num_trace
% [I,J,V]=find(Image(:,k));
% Image(:,k)=interp1q(J,V,1:1:NbPixZ);
%     % F=griddedInterpolant(I',Image(I,k),'nearest');
%     % Image(:,k)=F(X);
%     %Image(:,k)=griddata(I,J,V,Y,X);
%     %Image(:,k)=interp1((1:1:TimeWindow)',trace(77:end,k),(1:1:NbPixZ)','linear',0);
% end
% [I,J,V]=find(Image);
% [X,Y]=meshgrid(1:NbPixZ,1:NbPixX);
% Image=griddata(I,J,V,X,Y);
figure
imagesc(Image2);

%processing attenuation

%Mask outside
imagesc(diff(Image2,1,1));
diffImage2_smth=smooth2a(diff(Image2),10,10);
imagesc(diffImage2_smth);
logic2=diffImage2_smth~=0;
Image22=Image2.*[double(logic2);zeros(1,size(Image2,2))];
imagesc(Image22);

%compare 4 energy-average traces 
figure
slice=Image22(:,325);
slice=slice(slice~=0);
plot(slice);
hold all
slice=Image22(:,900);
slice=slice(slice~=0);
plot(slice);
slice=Image22(:,1900);
slice=slice(slice~=0);
plot(slice);
slice=Image22(:,2400);
slice=slice(slice~=0);
plot(slice);

diffImage2_smth=smooth2a(diff(Image2),250,250);
slicesmth=diffImage2_smth(:,2150);
slicesmth=slicesmth(slicesmth~=0);
plot(slicesmth);
















