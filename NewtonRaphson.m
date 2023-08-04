function [U,Epl,Sigma,Err]=NewtonRaphson(x,y,n,Btrans,A,N,W,u,epl,Sigma,BC,de_ThTr,deTrip,Tin,Tref,Xin,Xref,LGP,Err)%solution
tol=100;
aa1=1:length(u);
aa1(BC)=[];
count=1
DU=zeros(length(u),1);
du=zeros(length(u),1);
SigmaIn=Sigma;
Epl=epl;
Ftot=zeros(length(u),1);
while (tol>1e-3)
if count==1    
 [K , F, P, Epl, Sigma]=AssembleK(n,Btrans,A,W,u,du,DU,Epl,epl,Sigma,SigmaIn,de_ThTr,deTrip,Tin,Tref,Xin,Xref,LGP);
 Err=P+Err;
else
    Ftot=Ftot+F;
    Err=P-F;
%[K , F, P, Epl, Sigma]=AssembleK(n,Btrans,A,N,NT,W,u,du,epl,Sigma,BC,de_ThTr, deTrip,Tin,Tref,Xin,Xref,Tinitial1,Xinitial);
end
%K(length(u)*(BC-1)+BC)=1e10;
K(BC,:)=[]; %fixed boundary conditions applied on BC nodes
K(:,BC)=[];
Errsol=Err;
Errsol(BC)=[];
duPre=du;
[du1,flag,relres]=symmlq(K,Errsol,1e-6,1e4);
echo off
%du1=(K)\Errsol;
du=sparse(aa1,1,du1,length(u),1);
DU=DU+du;
U=u+DU;


%[K , F, P, Epl, Sigma]=AssembleK(n,Btrans,A,N,NT,W,u,DU,epl,Sigma,BC,de_ThTr, deTrip,Tin,Tref,Xin,Xref,Tinitial1,Xinitial);
[K , F, P, Epl, Sigma]=AssembleK(n,Btrans,A,W,u,du,DU,Epl,epl,Sigma,SigmaIn,de_ThTr,deTrip,Tin,Tref,Xin,Xref,LGP);
ErrPre=Err;
    Err=P-F;
ERR=Err;
ERR(BC)=[];
p=P;
p(BC)=[];
scatter3(x,y,F(1:2:end))
drawnow
hold on
scatter3(x,y,P(1:2:end))
drawnow
hold off
legend
%tol=((sqrt(sum(ERR.^2)))/sqrt(sum(p.^2)))
tol=abs((Err'*DU)./(P'*DU))
count=count+1
if count>200
tol=1e-4;
end
%if count>20

%end
end
