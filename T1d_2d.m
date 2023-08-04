function [Tstress_strain]=T1d_2d(T,x)
dx=max(abs(x(1:end-1)-x(2:end)));
xl=max(x);
xT=sort(unique(x));
[row col]=size(T);
 for i=1:length(xT)
    ii=find(x==xT(i));
     Tstress_strain(ii,1:col)=repmat(T(i,1:col),length(ii),1);
 end