function [E]=Boundry(x,y,x1,y1,n)
for count=1:length(x1)
[i(count),j]=find(and(x==x1(count),y==y1(count)));
end
[ii,jj]=size(n);
E=zeros(ii,jj+1);
for k=1:length(i)
[p q]=find(n==i(k));
for k2=1:length(p)
E(p(k2,1),1)=p(k2,1);
E(p(k2,1),q(k2,1)+1)=n(p(k2,1),q(k2,1));
end
end
[i j]=find(E(:,1)~=0);
E=E(i,:);
