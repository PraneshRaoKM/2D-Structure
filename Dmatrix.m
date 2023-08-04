function [D]= Dmatrix(E,mu)
D=E/((1+mu)*(1-2*mu))*[1-mu, mu, mu, 0; mu, 1-mu, mu, 0; mu, mu, 1-mu, 0; 0, 0, 0, (1-2*mu)/2 ];