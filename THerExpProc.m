
for ii=2:201
[row column]=size(n);
for i=1:row
    nodes=(n(i,:));
   v(i,1:2:8)=nodes*2-1;
   v(i,2:2:8)=nodes*2;

   FTherm(i,:)=ThermExp(Btrans(nodes,:,:),A(nodes,:),NT(nodes,:,:),T2d(nodes,ii-1));
   
end
P=sparse(reshape(v,[],1),1,reshape(FTherm,[],1),length(x)*2,1);
[u(:,ii)]=NewtonRaphson(x,y,n,u(:,ii-1),[BCnodes(E1,1,0);BCnodes(E2,0,1);BCnodes(E4,0,1)],P,T2d(:,ii),T2d(:,ii-1));
plot(u')

end