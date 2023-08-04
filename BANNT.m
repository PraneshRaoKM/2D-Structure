function [Btrans, A, N, NT, W]=BANNT(x,y,n,GP,GW)
r1=x(n(:,1));
r2=x(n(:,2));
r3=x(n(:,3));
r4=x(n(:,4));
z1=y(n(:,2));
z2=y(n(:,1));
z3=y(n(:,4));
z4=y(n(:,3));

[row ,~]=size(n);
count=1;
for i=1:length(GP)
    for j=1:length(GP)
Btrans(:,:,count)=BT(r1,r2,r3,r4,GP(i),GP(j),z1,z2,z3,z4);
A(:,count)=Pi2RDetJ(r1,r2,r3,r4,GP(i),GP(j),z1,z2,z3,z4);
N(:,:,count)=repmat(ShapeFunction(GP(i),GP(j)),row,1);
NT(:,:,count)=repmat(shapeFunT(GP(i),GP(j)),row,1);
W(count)=GW(i)*GW(j);
count=count+1;
    end
end
