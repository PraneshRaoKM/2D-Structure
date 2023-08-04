function [N]=shapeFunT(x_loc,y_loc)
N=[(x_loc-1.0).*(y_loc-1.0).*(1.0./4.0),(x_loc-1.0).*(y_loc+1.0).*(-1.0./4.0),(x_loc+1.0).*(y_loc+1.0).*(1.0./4.0),(x_loc+1.0).*(y_loc-1.0).*(-1.0./4.0)];