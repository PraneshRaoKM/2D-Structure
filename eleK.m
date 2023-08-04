function [K, F, P, S, eps]=eleK(BTrans,W,A,u,du,DU,T1,T2,X1,X2,de_Therm,deTripin,Sin,SigmaIn,Epl,epsin,LGP)
K=zeros(8,8);
F=zeros(8,1);
P=zeros(8,1);
for i=1:LGP
        BT=reshape(BTrans(1,:,i),[4,8]);
        %[D,dF,eps(:,:,i),S(:,:,i),Da]=NonLinearPlasticity(BT,T1(i),T1(i)-T2(i),X1(i,:)-X2(i,:),du,DU,X1(:,:,i),Sin(:,:,i),SigmaIn(:,:,i),epsin(:,:,i),de_Therm(i),deTripin(:,i));
         [D,dF,eps(:,:,i),S(:,:,i),Da]=ImplicitPlasticity(BT,T1(i),T1(i)-T2(i),X1(i,:)-X2(i,:),du,DU,X1(:,:,i),Sin(:,:,i),SigmaIn(:,:,i),Epl(:,:,i),epsin(:,:,i),de_Therm(i),deTripin(:,i));
       SS=Sin(:,:,i);
        Sdev=SS-[1,1,1,0]'*mean(SS(1:3));
        deTrip=Sdev*deTripin(:,i)';
        %deTrip=Sdev*0;
        de_ThTr=de_Therm(i)*[1 1 1 0]';
        dP=BT'*((D*(deTrip+de_ThTr)+Da')-SigmaIn);
        K= K+W(i)*A(i)*BT'*D*BT;

        F=F+A(i)*W(i)*dF;
        P=P+A(i)*W(i)*dP;
end