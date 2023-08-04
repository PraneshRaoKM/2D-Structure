for ii=240:nt
    ii
    Temp=T2d(:,ii);    
Tin(:,:,1:LGP)=sum(NT(:,:,1:LGP).*repmat(Temp(n),1,1,LGP),2);
XF=XF2d(:,ii);
XA=XA2d(:,ii);
XB=XB2d(:,ii);
XP=XP2d(:,ii);
XM=XM2d(:,ii);
Xin(:,:,1:LGP)=[sum(NT(:,:,1:LGP).*repmat(XF(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XA(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XB(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XP(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XM(n),1,1,LGP),2)];

[de_ThTr, deTrip]=QuenchLoadBC(Tin,Tref,Xin,Xref,Tinitial1,Xinitial,LGP);
BC=[BCnodes(E1,1,0);BCnodes(E2,0,1);BCnodes(E4,0,1)]; %Fixed BC
%[u(:,i),Epl(i,:,:,:),Sigma(i,:,:,:)]=NewtonRaphson(n,Btrans,A,N,W,u,epl,Sigma,BC,de_ThTr,deTrip,Tin,Tref,Xin,Xref)
[u(:,ii),epl(:,:,:,ii),Sigma(:,:,:,ii)]=NewtonRaphson(x,y,n,Btrans,A,N,W,u(:,ii-1),epl(:,:,:,ii-1),Sigma(:,:,:,ii-1),BC,de_ThTr, deTrip,Tin,Tref,Xin,Xref,LGP);%solution
subplot(1,2,1)
scatter3(x,y,u(1:2:end,ii))
subplot(1,2,2)
scatter3(x,y,u(2:2:end,ii))
drawnow
Tref=Tin;
Xref=Xin;
end