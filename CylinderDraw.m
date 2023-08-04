function [x, y, n, E1, E2, E3, E4]=CylinderDraw(xl,yl,dx,dy)
div_x=xl/dx; %number of divisions in in x direction
div_y=yl/dy; %number of divisions in in y direction
nx=div_x+1;%number of nodes in in x direction
ny=div_y+1;%number of nodes in in y direction
et=div_x*div_y; %Total numder of elements
nt=nx*ny;%Total numder of elements
%Node points
x=0:dx:xl;
y=0:dy:yl;
[X, Y]=meshgrid(x,y); % repeats x and y values to form grid
x=round(reshape(X,nt,1)/1000,6); %m to mm
y=round(reshape(Y,nt,1)/1000,6);
B=reshape(1:nt,ny,nx); %[1:nt] node number vector to grid
count=0;
n=zeros(et,4);

for j=1:div_y
    for i=1:div_x
    count=count+1;
    n(count,1:4)=reshape(B(j:j+1,i:i+1),1,4); %node connectivity matrix formation
    end
end
k=n;
n(:,1)=k(:,1);
n(:,2)=k(:,2);
n(:,3)=k(:,4);
n(:,4)=k(:,3);

E1=Boundry(x,y,x(find(x==0)),y(find(x==0)),n);
E2=Boundry(x,y,x(find(y==max(y))),y(find(y==max(y))),n);
E3=Boundry(x,y,x(find(x==max(x))),y(find(x==max(x))),n);
E4=Boundry(x,y,x(find(y==0)),y(find(y==0)),n);

clearvars -except x y n E1 E2 E3 E4 ntot dt theta beta epg;