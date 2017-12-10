%déconvolution

sample=reshape(stp(i:i+windowSize(1)-1,j:j+windowSize(2)-1),windowSize(1)*windowSize(2),1);
hPI = PI_deconvUknownth4(sample,'norm',Pnprimsigma(i,j)^2,Pnprimsigma(i,j));
fXdecUK=fdecUknown(xmesh,sample,10e10,'norm',Pnprimsigma(i,j),diff(xmesh(1:2)));
