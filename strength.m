function [H]= strength(Str,T,epsEq,TGP,HFO,HAO,HBO,HPO,HMO)
H=[interp2(Str,T,HFO',epsEq,TGP,'spline'),interp2(Str,T,HAO',epsEq,TGP,'spline'),interp2(Str,T,HBO',epsEq,TGP,'spline'),interp2(Str,T,HPO',epsEq,TGP,'spline'),interp2(Str,T,HMO',epsEq,TGP,'spline')];

  