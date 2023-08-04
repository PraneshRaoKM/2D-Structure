function [FTherm]=ThermExpBC(Btrans,A,NT,W,T,Tref,Tinitial)

FTherm=zeros(8,1);
for i=1:9
        BT=reshape(Btrans(:,:,i),[4,8]);
        TGP=NT(:,:,i)*T;
        TrefGP=NT(:,:,i)*Tref;
        D=Dmatrix(YM(TGP),PR(TGP));


        alphamat=[1 1 1 0]'*alpha(TGP)+[1 1 1 0]'*(alpha(TGP)-alpha(TrefGP))*(TGP-Tinitial)/(TGP-TrefGP);
        dF=alphamat'*A(i)*D*BT*(TGP-TrefGP);
        FTherm=FTherm+W(i)*dF';
end