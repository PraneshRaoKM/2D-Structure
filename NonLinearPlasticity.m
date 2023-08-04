function [D,Force,eps,S,Da]=NonLinearPlasticity(BT,TGP,dT,dX,du,DU,XGP,Sin,SigmaIn,epsin,deThin,deTripin)
load('HMechProp.mat');
%YEIELD-STRENGTH and H
de=BT*DU;
count=1;
Ephase=interp1(Tprop,E,TGP);
YM=Ephase*XGP';
PoiRat=interp1(Tprop,PR,TGP)*XGP';
Del=Dmatrix(YM,PoiRat);
Da=zeros(1,4);
Sdev=Sin-[1,1,1,0]'*mean(Sin);
%deTrip=Sdev*deTripin;
deTrip=Sdev*deTripin;
de_ThTr=deThin*[1 1 1 0]';
Strial=SigmaIn+Del*(de-deTrip-de_ThTr);


epsEq=sqrt(2/3*epsin'*epsin);
YS_Phase=interp1(T,YS,TGP);
YS_T=YS_Phase*XGP';
YS_dT=((YS_T-interp1(T,YS,.99*TGP))*XGP')/(.01*TGP);
YS_dX=YS_Phase*dX';

H_phase=strength(Str,T,epsEq,TGP,HFO,HAO,HBO,HPO,HMO);
%H_phase=((strength(YS_Phase,Ephase,Strain,T,1.01*(YS_T/YM),TGP,SFer,SAus,SBain,SPaer,SMar))-strength(Strain,T,(YS_T/YM),TGP,SFer,SAus,SBain,SPaer,SMar))/(.01*(YS_T/YM));
H_T=H_phase*XGP';
H_dT=((H_T-strength(Str,T,epsEq,TGP,HFO,HAO,HBO,HPO,HMO))*XGP')/(.01*TGP);
H_dX=H_phase*dX';
H_depsq=((H_T-strength(Str,T,.99*epsEq,TGP,HFO,HAO,HBO,HPO,HMO))*XGP')/(.01*epsEq);

[Seq, df_ds, d2_f2]=equivalent(Strial);
%Seq=sqrt(((Strial(1)-Strial(2))^2+(Strial(2)-Strial(3))^2+(Strial(3)-Strial(1))^2+6*Strial(4)^2)/2);
f=Seq-(YS_T+H_T);
if f<0
    D=Del;
    eps=epsin;
    S=Strial;
    
else

        S=Strial;
    DY=0;
    depsEq=0;
    
   %df_ds=(1/(4*Seq))*[(4*S(1)-2*S(2)-2*S(3));(4*S(2)-2*S(1)-2*S(3));(4*S(3)-2*S(1)-2*S(2));12*S(4)];
    Deps=zeros(size(epsin));
while or(abs(f)>1e3,count==1) 
YS_Phase=interp1(T,YS,TGP);
YS_T=YS_Phase*XGP';
YS_dT=((YS_T-interp1(T,YS,.99*TGP))*XGP')*dT/(.01*TGP);
YS_dX=YS_Phase*dX';

H_phase=strength(Str,T,epsEq,TGP,HFO,HAO,HBO,HPO,HMO);
%H_phase=((strength(YS_Phase,Ephase,Strain,T,1.01*(YS_T/YM),TGP,SFer,SAus,SBain,SPaer,SMar))-strength(Strain,T,(YS_T/YM),TGP,SFer,SAus,SBain,SPaer,SMar))/(.01*(YS_T/YM));
H_T=H_phase*XGP';
H_dT=((H_T-strength(Str,T,1.01*epsEq,TGP,HFO,HAO,HBO,HPO,HMO))*XGP')*dT/(.01*TGP);
H_dX=H_phase*dX';
if epsEq==0
    epsEq=YS_T/YM;
end
    H_depsq=((H_T-strength(Str,T,1.01*epsEq,TGP,HFO,HAO,HBO,HPO,HMO))*XGP')/(.01*epsEq);
    
    dY=f/((df_ds'*Del*df_ds)); %change values of H with updated eps    
    if dY<0
  dY=0;
    end
    DY=DY+.5*dY;
    S=S-.5*dY*Del*df_ds;
    [Seq, df_ds, d2_f2]=equivalent(S);
    %Seq=sqrt(((S(1)-S(2))^2+(S(2)-S(3))^2+(S(3)-S(1))^2+6*S(4)^2)/2);
    %df_ds=(1/(4*Seq))*[(4*S(1)-2*S(2)-2*S(3));(4*S(2)-2*S(1)-2*S(3));(4*S(3)-2*S(1)-2*S(2));12*S(4)];
    deps=dY*df_ds;
    Deps=Deps+deps;
    eps=epsin+Deps;
    epsEq=sqrt(2/3*eps'*eps); 
    
    H_phase=strength(Str,T,epsEq,TGP,HFO,HAO,HBO,HPO,HMO);
    H_T=H_phase*XGP';
    f=(Seq)-(YS_T+H_T); %change values of YS with updated eps
    count=count+1;
     xxx(count)=f;
     yyy(count)=YS_T+H_T*epsEq;
     %plot(xxx)
     %drawnow
     if count==15
         break
     end
end
    sqrt(2/3*df_ds'*df_ds)
    D=Del-(Del*(df_ds'*Del*df_ds)/((df_ds'*Del*df_ds)+sqrt(2/3*df_ds'*df_ds)*H_depsq));
    Da=(-df_ds'*Del*(H_dT+H_dX))/((df_ds'*Del*df_ds)+sqrt(2/3*df_ds'*df_ds)*H_depsq);
    %D=Del-Del*(df_ds'*Del*df_ds)/(H_depsq+(df_ds'*Del*df_ds)); %change values of H with updated eps
    %Sdev=S-[1,1,1,0]'*mean(S(1:3));
    %deTrip=Sdev*deTripin*1e-6;
end
%Force=BT'*D*(de-deTrip-de_ThTr);
Force=BT'*(S-SigmaIn);


%Force=BT'*(D*(de-deTrip-de_ThTr)+Da');
