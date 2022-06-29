function [numflux] = LinwaveLF(u,v,lambda,wavesp)
% function [numflux] = LinwaveLF(u,v,lambda,maxvel);
% Purpose: Evaluate Lax Friedrich numerical flux for wave equation

numflux = wavesp*(u+v)/2 - abs(wavesp)/2*(v-u);
end