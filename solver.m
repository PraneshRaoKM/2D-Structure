clear
load('D15_h1500_Tb350.mat')
[x, y, n, E1, E2, E3, E4]=CylinderDraw(7.5,2,.25,.5); %2D cylinder mesh generation
GP=[0];
GW=[2];
%GP=[sqrt(3/5) 0 -sqrt(3/5)];
%GW=[5/9 8/9 5/9];
[~, nt]=size(T);
[T2d]=T1d_2d(interp1((0:.1:(length(T)-1)*.1),T',(0:.1:(length(T)-1)*.1))',x); %Temperature profile from 1D to 2D
[XF2d]=T1d_2d(interp1((0:.1:(length(T)-1)*.1),XF',(0:.1:(length(T)-1)*.1))',x);
[XB2d]=T1d_2d(interp1((0:.1:(length(T)-1)*.1),XB',(0:.1:(length(T)-1)*.1))',x);
[XP2d]=T1d_2d(interp1((0:.1:(length(T)-1)*.1),XP',(0:.1:(length(T)-1)*.1))',x);
[XM2d]=T1d_2d(interp1((0:.1:(length(T)-1)*.1),XM',(0:.1:(length(T)-1)*.1))',x);
[XA2d]=1-(XF2d+XB2d+XP2d+XM2d);

[Btrans, A, N, NT, W]=BANNT(x,y,n,GP,GW);%node and shape matrices
LGP=length(GP)*length(GP);
Temp=T2d(:,1);   
Tinitial1(:,:,1:LGP)=sum(NT(:,:,1:LGP).*repmat(Temp(n),1,1,LGP),2);
Tref=Tinitial1;
%Tinitial1(:,:,1:LGP)=repmat(25,size(Tinitial1));
%Tref=Tinitial1;
XF=XF2d(:,1);
XA=XA2d(:,1);
XB=XB2d(:,1);
XP=XP2d(:,1);
XM=XM2d(:,1);
Xinitial(:,:,1:LGP)=[sum(NT(:,:,1:LGP).*repmat(XF(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XA(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XB(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XP(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XM(n),1,1,LGP),2)];
Xref=[sum(NT(:,:,1:LGP).*repmat(XF(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XA(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XB(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XP(n),1,1,LGP),2),sum(NT(:,:,1:LGP).*repmat(XM(n),1,1,LGP),2)];
Xref=Xinitial;
Xinitial(:,:,1:LGP)=repmat([0 0 0.75 .25 0],120,1);

u=zeros(length(x)*2,1); %initial condition
epl=zeros(4,length(n),LGP,1);
Sigma=zeros(4,length(n),LGP,1);
Err=zeros(length(x)*2,1);
for ii=2:nt*10
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
[u(:,ii),epl(:,:,:,ii),Sigma(:,:,:,ii), Err]=NewtonRaphson(x,y,n,Btrans,A,N,W,u(:,ii-1),epl(:,:,:,ii-1),Sigma(:,:,:,ii-1),BC,de_ThTr, deTrip,Tin,Tref,Xin,Xref,LGP,Err);%solution
subplot(1,2,1)
scatter3(x,y,u(1:2:end,ii))
subplot(1,2,2)
scatter3(x,y,u(2:2:end,ii))
drawnow
Tref=Tin;
Xref=Xin;
end