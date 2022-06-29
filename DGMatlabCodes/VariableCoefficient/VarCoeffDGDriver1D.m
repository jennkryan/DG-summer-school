% Driver script for solving the 1D Variable Coefficient equation using an 
% DG scheme
clear all
close all

% Order of method (m), number of elements (N)
m=16; 
N=5;
    
% Set problem parameters
xmin = -1.0; xmax = 1.0; L=xmax - xmin;
FinalTime = 10.5;
CFL = 0.1;

% Generate mesh
VX = (xmax-xmin)*(0:N)/N + xmin; 
r = LegendreGL(m);
x = ones(m+1,1)*VX(1:N) + 0.5*(r+1)*(VX(2:N+1)-VX(1:N));
h = (xmax-xmin)/N;

% Define initial conditions
time = 10.5;
u = sin(4*pi*x); %periodic BC needed
% u = (1-sign(x-0.2))/2+1; % Constant BC needed

% Define variable coefficient
Avar = (1-x.^2).^5+1;

% Solve Problem
[u] = VarCoeffDG1D(x,Avar,u,h,m,N,CFL,FinalTime);

plot(x,u,'-o','LineWidth',3)
xlabel('X','FontSize', 14)
ylabel('u_h(x,T)','FontSize', 14)
title(['Approximation for p=',num2str(m),...
    ' and N=',num2str(N),' elements'],'FontSize', 16)