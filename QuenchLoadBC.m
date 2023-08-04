function [de_ThTr, deTrip]=QuenchLoadBC(Tin,Tref,Xin,Xref,Tinitial,Xinitial,LGP)
de_ThTr=zeros(length(Tin),LGP);
deTrip=zeros(length(Tin),LGP);
for i=1:LGP     
    TGP=Tin(:,:,i);
    TrefGP=Tref(:,:,i);
    TinitialGP=Tinitial(:,:,i);
    XGP=Xin(:,:,i);
    XrefGP=Xref(:,:,i);
    XinitialGP=Xinitial(:,:,i);
    
    load('Mechprop.mat');
        
        rhoRef=sum(interp1(Tprop,Rho,TinitialGP).*XinitialGP,2);
        rho1=sum(interp1(Tprop,Rho,TGP).*XGP,2);
        rho2=sum(interp1(Tprop,Rho,1.01*TGP).*XGP,2);        
        drho_T=((rho2-rho1))./(.01*TGP);
        dT=(TGP-TrefGP);
        drho_t=drho_T.*dT;
        dX=((XrefGP-XGP));
        drho_X=sum(interp1(Tprop,Rho,TGP).*dX,2);  

        de_dRho=-1/3*((rhoRef).^(1/3)).*(rho1.^(-4/3));

       de_ThTr(:,i)=de_dRho.*(drho_t+drho_X);
        
      
        XGP(:,2)=[];
        dX(:,2)=[];
        [ss,vv]=find(XGP>0);
        GW=interp1(Tprop,1e-6*TransPlas,TGP);
        deTrip(ss,i)=deTrip(ss,i)-3/2*sum(GW(ss,vv).*dX(ss,vv).*log(XGP(ss,vv)),2);

end

        