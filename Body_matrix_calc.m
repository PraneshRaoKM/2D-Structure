N1=@(x_loc,y_loc)(1-x_loc)*(1-y_loc)/4;
N2=@(x_loc,y_loc)(1-x_loc)*(1+y_loc)/4;
N3=@(x_loc,y_loc)(1+x_loc)*(1+y_loc)/4;
N4=@(x_loc,y_loc)(1+x_loc)*(1-y_loc)/4;
syms x_loc y_loc r1 r2 r3 r4 z1 z2 z3 z4 dx dy
assume([x_loc y_loc r1 r2 r3 r4 z1 z2 z3 z4 dx dy],'real')
N=[N1(x_loc,y_loc) 0 N2(x_loc,y_loc) 0 N3(x_loc,y_loc) 0 N4(x_loc,y_loc) 0;0 N1(x_loc,y_loc) 0 N2(x_loc,y_loc) 0 N3(x_loc,y_loc) 0 N4(x_loc,y_loc)];
N_J=[N1(x_loc,y_loc) N2(x_loc,y_loc) N3(x_loc,y_loc) N4(x_loc,y_loc)];
r=N_J*[r1 r2 r3 r4]';
z=N_J*[z1 z2 z3 z4]';
J=[diff(r,x_loc),diff(z,x_loc);diff(r,y_loc),diff(z,y_loc)];
invJ=inv(J);
INVJ=[invJ(1,1),invJ(2,1),0, 0; invJ(2,1),invJ(2,2),0, 0; 0, 0, invJ(1,1),invJ(2,1); 0, 0, invJ(2,1),invJ(2,2)];
A1=[1 0 0 0; 0 0 0 1; 0 0 0 0; 0 1 1 0]*INVJ;
G=[diff(N(1,:),x_loc); diff(N(1,:),y_loc); diff(N(2,:),x_loc);diff(N(2,:),y_loc)];
A2=[0, 0; 0, 0; 1/r, 0; 0, 0;];
B=(A1*G+A2*N);
funG=matlabFunction(simplify(G));
funB=matlabFunction(simplify(B));
funN=matlabFunction(simplify(N));