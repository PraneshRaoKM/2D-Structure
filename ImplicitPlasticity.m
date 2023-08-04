function [D,Force,eps,S,Da]=ImplicitPlasticity(BT,TGP,dT,dX,du,DU,XGP,Sin,SigmaIn,Epl,epsin,deThin,deTripin)
load('HMechProp.mat');
%YEIELD-STRENGTH and H
eps=epsin;
de=BT*DU;
epsinEq=sqrt(2/3*epsin'*epsin);
count=1;
Ephase=interp1(Tprop,E,TGP);
YM=Ephase*XGP';
PoiRat=interp1(Tprop,PR,TGP)*XGP';
Del=Dmatrix(YM,PoiRat);
Da=zeros(1,4);

de_ThTr=deThin*[1 1 1 0]';
II=[2/3 -1/3 -1/3 0; -1/3 2/3 -1/3 0;-1/3 -1/3 2/3 0; 0 0 0 1];
Strial=inv(eye(4)+Del*deTripin*II)*(SigmaIn+Del*(de-de_ThTr));
 

epsEq=epsinEq;
YS_Phase=interp1(T,YS,TGP);
YS_dT=-(YS_Phase-interp1(T,YS,1.01*TGP))*XGP'/(.01*TGP);
YS_T=YS_Phase*XGP';

H_phase=strength(Str,T,epsEq,TGP,HFO,HAO,HBO,HPO,HMO);
H_T=H_phase*XGP';
H_dT=((H_T-strength(Str,T,epsEq,TGP,HFO,HAO,HBO,HPO,HMO))*XGP')/(.01*TGP);
H_depsq=-((H_T-strength(Str,T,epsEq+1e-6,TGP,HFO,HAO,HBO,HPO,HMO))*XGP')/(1e-6);
[Seq, df_ds, d2_f2]=equivalent(Strial);
f=Seq-(YS_T+H_T);
if f<5000
    D=Del;
    eps=epsin;
    S=Strial;
    if Strial>0
        w=(YS_T+H_T)/Seq;
        S=w*S;
    end
        
else

        S=Strial;
    DY=0;
    depsEq=0;
    
   %df_ds=(1/(4*Seq))*[(4*S(1)-2*S(2)-2*S(3));(4*S(2)-2*S(1)-2*S(3));(4*S(3)-2*S(1)-2*S(2));12*S(4)];
    Deps=zeros(size(epsin));


dY=0;
f=Seq-(YS_T+H_T);
eps=epsin;
d_Sdev=[2/3 -1/3 -1/3 0;-1/3 2/3 -1/3 0; -1/3  -1/3 2/3 0;0 0 0 1];
I=eye(4);
tol=100;
counter=1;
P=zeros(4,1);
counter=1;
while (f)>5000
    toler=tol;
d_P_ds=I+Del*deTripin*d_Sdev+Del*dY*d2_f2;
%d_P_ds=I+Del*dY*d2_f2;
d_P_dY=Del*df_ds;
dF_ds=df_ds';
dF_dY=H_depsq*(2/3)*sqrt(df_ds'*df_ds);
J=[d_P_ds d_P_dY; dF_ds dF_dY];
delta=J\[-P; -f];
if delta(5)<0
    break;
end
dSig=.5*delta(1:4);
S=S+dSig;
[Seq, df_ds, d2_f2]=equivalent(S);
dY=dY+.5*delta(5);
P=(S-SigmaIn)-Del*(de-de_ThTr-deTripin*d_Sdev*S)+Del*dY*df_ds;
[dSeq, ~, ~]=equivalent(delta(1:4));

deps=dY*df_ds;
eps=epsin+deps;
epsEq=sqrt(2/3*eps'*eps);

depsEq=sqrt(2/3*deps'*deps);
H_phase=strength(Str,T,epsEq,TGP,HFO,HAO,HBO,HPO,HMO);
H_T=H_phase*XGP';

H_depsq=-((H_T-strength(Str,T,epsEq+1e-6,TGP,HFO,HAO,HBO,HPO,HMO))*XGP')/(1e-6);
f=Seq-(YS_T+H_T);
%tol=max(abs([(dSeq/Seq),(dY/delta(5))]));
counter=counter+1;
end
w=(YS_T+H_T)/Seq;
S=w*S;
[Seq1, df_ds1, d2_f21]=equivalent(S);
[Seq, df_ds, d2_f2]=equivalent(SigmaIn+.0001);

H_dX=H_phase*dX';
YS_dX=YS_Phase*dX';
Den=(df_ds1'*Del*df_ds1)+(2/3*sqrt(df_ds1'*df_ds1)*H_depsq);
 D=Del-(Del*(df_ds1*(df_ds1'*Del)))/Den;
 Da=(df_ds1'*Del*(YS_dX+YS_dT++H_dT+H_dX))/Den;
  %Sdev=S-[1,1,1,0]'*mean(S(1:3));
 %deTrip=Sdev*deTripin;
 %dF=(D*(de-deTrip-de_ThTr)+Da');
%Force=((S-SigmaIn));
%dF./Force;
end
 %Sdev=S-[1,1,1,0]'*mean(S(1:3));
 %deTrip=Sdev*deTripin;
 %dS=(D*(de-deTrip-de_ThTr)+Da');
 %S=SigmaIn+dS;
%Force=BT'*dS;
Force=BT'*((S-SigmaIn));