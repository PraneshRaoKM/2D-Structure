function [de_ThTr, deTrip]=ThermExpBC(n,Btrans,A,NT,W,Temp,Tref,X,Xref,Tinitial,Xinitial)
[row column]=size(n);
    
for i=1:row
    nodes=n(i,:);
    
  [de_ThTr(i,:), deTrip(i,:)]=QuenchLoadBC(Btrans(i,:,:),A(i,:),NT(i,:,:),W,Temp(nodes,:),Tref(nodes,:),X(nodes,:),Xref(nodes,:),Tinitial(nodes,:),Xinitial(nodes,:));
end