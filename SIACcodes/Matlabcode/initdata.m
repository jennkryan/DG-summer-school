

evalPoints = input('Input number of evaluation points per element:  ');

% this WILL break if you don't have enough elements to support the degree
% of the kernel.  If you receive an index out of bounds error, increase the
% number of elements
N=input('Input number of elements:  ');  %number of elements
M = input('Input polynomial degree (NOT ORDER):  ');   %polynomial degree

%%%%% UNCOMMENT THESE TWO LINES FOR UNIFORM INTERVALS
xlen = xmax-xmin;
h = xlen/N;
x = xmin:h:xmax;

%%%%% UNCOMMENT THIS SINGLE LINE FOR NON_UNIFORM INTERVAL
% x = nonUniformIntervals( N, 1.0, 0 ); %element definitions

[z,w] = JacobiGZW(M+1,0.0,0.0);  %get quadrature points and weights
% six quadrature points at which to evaluate the funcitons -- non
% superconvergent
[zEval, wEval] = JacobiGZW(evalPoints, 0.0, 0.0);
NumP = size(zEval, 1);


Lvalue = zeros(M+1);
for i = 0:M
  Lvalue(i+1,:) = sqrt(i+0.5)*JacobiPoly(i,z,0.0,0.0)';
  LvalueMass(i+1) = dot(w,Lvalue(i+1,:).*Lvalue(i+1,:));
end



zmapquad = ones(M+1,1)*x(1:N) + 0.5*(z+1)*(x(2:N+1)-x(1:N));
zmapEval = ones(NumP,1)*x(1:N) + 0.5*(zEval+1)*(x(2:N+1)-x(1:N));
