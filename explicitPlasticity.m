function [D,Force,eps,S,Da]=explicitPlasticity(BT,TGP,dT,dX,du,DU,XGP,Sin,SigmaIn,epsin,deThin,deTripin)
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
if f<0
    D=Del;
    eps=epsin;
    S=Strial;        
else

S=Strial;
f=Seq-(YS_T+H_T);
eps=epsin;
d_Sdev=[2/3 -1/3 -1/3 0;-1/3 2/3 -1/3 0; -1/3  -1/3 2/3 0;0 0 0 1];
I=eye(4);
tol=100;
counter=1;
P=zeros(4,1);
counter=1;
w=(YS_T+H_T)/Seq;
S=w*S;
[Seq, ~, d2_f2]=equivalent(S);
[dSeq, df_ds, d2_f2]=equivalent(Strial-SigmaIn);
[SeqIn, df_ds1, d2_f21]=equivalent(SigmaIn);

H_dX=H_phase*dX';
YS_dX=YS_Phase*dX';
Den=(df_ds1'*Del*df_ds1)+(2/3*sqrt(df_ds1'*df_ds1)*H_depsq);
 D=Del-(Del*(df_ds1*(df_ds1'*Del)))/Den;
 Da=(df_ds1'*Del*(YS_dX+YS_dT++H_dT+H_dX))/Den;
  Sdev=S-[1,1,1,0]'*mean(S(1:3));
 deTrip=Sdev*deTripin;
 dF=(D*(de-deTrip-de_ThTr)+Da');
Force=((S-SigmaIn));
dF./Force
end
 %Sdev=S-[1,1,1,0]'*mean(S(1:3));
 %deTrip=Sdev*deTripin;
 %dS=(D*(de-deTrip-de_ThTr)+Da');
 %S=SigmaIn+dS;
%Force=BT'*dS;
Force=BT'*((S-SigmaIn));


