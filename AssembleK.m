function [K , F, P, eps, S]=AssembleK(n,Btrans,A,W,u,du,DU,Epl,epl,Sigma,SigmaIn,de_ThTr, deTrip,T1,T2,X1,X2,LGP)
l=length(u); %2 DOF per node
[row, ~]=size(n);
for i=1:row
Knodes(1:2:8,1)=2*n(i,:)-1;
Knodes(2:2:8,1)=2*n(i,:);

[Kele, Fele, Pele, S(:,i,:), eps(:,i,:)]=eleK(Btrans(i,:,:),W,A(i,:),u(Knodes),du(Knodes),DU(Knodes),T1(i,1,:),T2(i,1,:),X1(i,:,:),X2(i,:,:),de_ThTr(i,:),deTrip(i,:),Sigma(:,i,:),SigmaIn(:,i,:),Epl(:,i,:),epl(:,i,:),LGP);
%[Kele, Fele]=eleK(Btrans(i,:,:),W,A(i,:),u(Knodes),du(Knodes),T(n(i,:)));
[nn1, nn2, v]=find(Kele);

kv1(i,:)=Knodes(nn1);
kv2(i,:)=Knodes(nn2);
vec(i,:)=v;
kvF(i,:)=Knodes;
vecF(i,:)=Fele;
vecP(i,:)=Pele;
end
K=sparse(reshape(kv1,[],1),reshape(kv2,[],1),reshape(vec,[],1),l,l);
F=sparse(reshape(kvF,[],1),1,reshape(vecF,[],1),l,1);
P=sparse(reshape(kvF,[],1),1,reshape(vecP,[],1),l,1);