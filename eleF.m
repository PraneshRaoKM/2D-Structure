function [F]=eleF(r,z,u,T,Tref)
GP=[sqrt(3/5) 0 -sqrt(3/5)];
GW=[5/9 8/9 5/9];
K=zeros(8,8);
eps=zeros(4,1);
sig=zeros(4,1);
F=zeros(8,1);
for i=1:3
    for j=1:3
        BTrans=BT(r(1),r(2),r(3),r(4),GP(i),GP(j),z(1),z(2),z(3),z(4));
        W=GW(i)*GW(j);
        A=Pi2RDetJ(r(1),r(2),r(3),r(4),GP(i),GP(j),z(1),z(2),z(3),z(4));
        TGP=shapeFunT(GP(i),GP(j))*T;
        D=Dmatrix(YM(TGP),PR(TGP));
        deps=BTrans*u;
        
        alphamat=[1 1 1 0]'*alpha(TGP);
        TrefGP=(shapeFunT(GP(i),GP(j)))*Tref;
        deth=alphamat'*(TGP-TrefGP);
        
        dSigma=D*(deps);
        dF=dSigma'*BTrans;
        eps=eps+W*deps;
        sig=sig+GW(i)*GW(j)*dSigma;
        F=F+A*W*dF';
    end
end