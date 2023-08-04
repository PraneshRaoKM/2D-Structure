function [eps sig F]=stress_strain(r,z,u,D)
GP=[sqrt(3/5) 0 -sqrt(3/5)];
GW=[5/9 8/9 5/9];
eps=zeros(4,1);
sig=zeros(4,1);
F=zeros(8,1);
for i=1:3
    for j=1:3
    deps=BT(r(1),r(2),r(3),r(4),GP(i),GP(j),z(1),z(2),z(3),z(4))*u
    dSigma=D*deps
    dF=dSigma'*BT(r(1),r(2),r(3),r(4),GP(i),GP(j),z(1),z(2),z(3),z(4));
    eps=eps+GW(i)*GW(j)*deps
    sig=sig+GW(i)*GW(j)*dSigma
    F=F+GW(i)*GW(j)*dF'
    end
end
