function [numflux] = VarCoeffSFCF(BvarU,BvarV,u,v,lambda,maxvel)
% function [numflux] = VarCoeffLF(BvarU,BvarV,u,v,lambda,maxvel);
% Purpose: Evaluate the Lax Friedrich numerical flux for Burgers equation
% XXXX USING CENTRAL FLUX?  I THINK THERE WAS A TYPO IN HERE.

%Burgers
%fu = 0.5.*u.^2; fv = 0.5.*v.^2;
% Linear case:  fu = u
%fu = u; fv = v;

%Variable Coefficient
fu = BvarU.*u; fv = BvarV.*v;

%numflux = 0.5*(fu-fv);

numflux = 0.5*(fu+fv);

end
