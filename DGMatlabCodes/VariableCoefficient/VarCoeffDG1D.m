function [u] = VarCoeffDG1D(x,Avar,u,h,m,N,CFL,FinalTime)
% function [u] = VarCoeffDG1D(x,Avar,u,h,m,N,CFL,FinalTime)
% Purpose  : Integrate 1D variable coefficient equation until FinalTime using a DG
%            scheme and 3rd order SSP-RK method
% Initialize operators at Legendre Gauss Lobatto grid
r = LegendreGL(m); 
V  = VandermondeDG(m, r); 
D = DmatrixDG(m, r, V);
Ma = inv(V*V'); 
S = Ma*D; 
iV = inv(V);

% Compute operator for WENO smoothness evaluator
%[Q,Xm,Xp] = WENODGWeights(m,iV);

% Initialize extraction vector
VtoE = zeros(2,N);
for j=1:N
  VtoE(1,j) = (j-1)*(m+1)+1; % first node in element j
  VtoE(2,j) = j*(m+1); % last node in element j
end

% Initialize filter matrix
% F = FilterDG(m,0,10,V);

% Compute smallest spatial scale timestep
rLGLmin = min(abs(r(1)-r(2))); 
time = 0; 
tstep = 0;

% Initialize parameters for nonlinear viscosity
% nu = zeros(m+1,N); 
% nu0 = 2; 
% kappa = -6; 
% c2 = 1;

maxvel = max(max(abs(Avar))); 

% integrate scheme
while (time<FinalTime)
  % Decide on timestep
  
  k = CFL*rLGLmin*h/maxvel; 
  if (time+k>FinalTime) 
      k = FinalTime-time; 
  end  
  % Update solution - stage 1
  rhsu  = VarCoeffDGrhs1D(x,Avar,u,h,k,m,N,Ma,S,VtoE,maxvel,time); 
  u1 = u + k*rhsu;
  %u1 = WENOlimitDG(x,u1,m,h,N,V,iV,Q,Xm,Xp);
  %u1 = MomentLimitDG(x,u1,m,h,N,V,iV);
  
  % Update solution - stage 2
  rhsu  = VarCoeffDGrhs1D(x,Avar,u1,h,k,m,N,Ma,S,VtoE,maxvel,time);
  u2 = (3*u + u1 + k*rhsu)/4;
  %u2 = WENOlimitDG(x,u2,m,h,N,V,iV,Q,Xm,Xp);
  %u2 = MomentLimitDG(x,u2,m,h,N,V,iV);
  
  % Update solution - stage 3
  rhsu  = VarCoeffDGrhs1D(x,Avar,u2,h,k,m,N,Ma,S,VtoE,maxvel,time); 
  u = (u + 2*u2 + 2*k*rhsu)/3;
  %u = WENOlimitDG(x,u,m,h,N,V,iV,Q,Xm,Xp);
  %u = MomentLimitDG(x,u,m,h,N,V,iV);
  time = time+k; 
  tstep = tstep+1;
end
return