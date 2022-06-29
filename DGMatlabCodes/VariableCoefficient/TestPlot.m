clear all
close all

% Order of method (m), number of elements (N)
m=16; 
N=5;

% Set problem parameters
xmin = -1.0; xmax = 1.0; L=xmax - xmin;

% Generate mesh
VX = (xmax-xmin)*(0:N)/N + xmin; 
r = LegendreGL(m);
x = ones(m+1,1)*VX(1:N) + 0.5*(r+1)*(VX(2:N+1)-VX(1:N));
h = (xmax-xmin)/N;

u = sin(4*pi*x);

Avar = (1-x.^2).^5+1;

Y = zeros((m+1)*N,1);
Aplot = Y;
for j=1:N
    for k=1:m+1
        Y((j-1)*(m+1)+k) = x(k,j);
        Aplot((j-1)*(m+1)+k) = Avar(k,j);
    end
end

plot(Y,Aplot)
