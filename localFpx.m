FpxCell=cell(1,length(v));

for k=1:length(v)
    k=7;
    %(5) Migrate this superposition over the previous range of velocities
    [s_hat_mig{k},tshatmig_n,xshatmig_n]=fkmig(s_hat,0.8e-9,0.2,v(k));
    
    N = size(s_hat_mig{k},1) - 2*windowSize(1) + 1;
    M = size(s_hat_mig{k},2) - 2*windowSize(2) + 1;
    FpxCell{k}=zeros(N,M);

    for i = 1:N
        for j = 1:M
            %estimate P(x) locally from kernel density estimation.
            sample=s_hat_mig{k}(i:i+windowSize(1)-1,j:j+windowSize(2)-1);
            %[~,Px,xmesh,~] = kde(sample,2^8,-10e4,10e4);
            [~,Px,xmesh,~] = kde(sample,2^8);

            
            Px=abs(Px');
            I=find(Px==0);
            while isempty(I)==0
                I=find(Px==0);
                for l=1:length(I)
                    if I(l)==1
                        Px(I(l))=Px(I(l)+1);
                        
                    elseif I(l)==length(Px)
                        Px(I(l))=Px(I(l)-1);
                        
                    else
                        Px(I(l))=0.5*(Px(I(l)-1)+Px(I(l)+1));
                        
                    end
                    
                end
                I=find(Px==0);
            end
            
            FpxCell{k}(i,j)=sum(Px.*log10(Px)*diff(xmesh(1:2)))+0.5*log10(sum((xmesh.^2).*Px*diff(xmesh(1:2))))+log10(sqrt(2*pi))+0.5;
        end
    end
end

CellFpx=cell(N,M);
for k=1:length(v)
for i=1:N
    for j=1:M
        CellFpx{i,j}(k)=FpxCell{k}(i,j);
    end
end
end

figure
for k=1:length(v)
subplot(3,4,k)
imagesc(FpxCell{k},[0 0.8])
end