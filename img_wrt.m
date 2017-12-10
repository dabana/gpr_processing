function [Image] = img_wrt(time,trace,ElevMat,NbPixZ,resoZ,NbPixX,resoX,ele0)
%Fonction qui écrit une image
%
Image=NaN(NbPixZ,NbPixX);
TimeWindow=length(time(66:end));
  for i=1:NbPixZ
      diff=abs(ElevMat-(ele0-i*resoZ*ones(TimeWindow,NbPixX)));
      logic=diff==repmat(min(diff),TimeWindow,1);
      Image(i,:)=sum(logic.*trace(66:length(time),:),1);
      clear diff logic
  end

end

