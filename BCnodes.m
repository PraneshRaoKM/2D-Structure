function [BC]=BCnodes(E,x,y)
nodes=unique(E(find(E(:,2:end)>0)+length(E(:,1))));
if and((x==1),(y==1))
    BC=[nodes*2-1;nodes*2];
elseif and((x==1),(y==0))
    BC=nodes*2-1;
elseif and((x==0),(y==1))
    BC=nodes*2;
end