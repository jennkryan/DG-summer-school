function [numflux] = VarCoeffLF(BvarU,BvarV,u,v,lambda,maxvel)
% function [numflux] = VarCoeffLF(BvarU,BvarV,u,v,lambda,maxvel);
% Purpose: Evaluate the Lax Friedrich numerical flux for Burgers equation

%Burgers
%fu = 0.5.*u.^2; fv = 0.5.*v.^2;
% Linearr case:  fu = u
%fu = u; fv = v;
%Variable Coefficient
fu = BvarU.*u; fv = BvarV.*v;

numflux = (fu+fv)/2 - maxvel/2*(v-u);

end