function [YSF]=yield(Strain,SF,E,Tprop)
a=[25 50:50:900];
e=interp1(Tprop,E,a);
m=(SF(2,:)-SF(1,:))/Strain(2);
size(m)
size(e)
YSStrF=SF(1,:)./(e-(m));
for i=1:19
YSF(i)=interp1(Strain(:,i),SF(:,i),YSStrF(i));
end