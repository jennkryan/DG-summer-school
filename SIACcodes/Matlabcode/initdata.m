

evalPoints = input('Input number of evaluation points per element:  ');

% this WILL break if you don't have enough elements to support the degree
% of the kernel.  If you receive an index out of bounds error, increase the
% number of elements
N=input('Input number of elements:  ');  %number of elements
M = input('Input polynomial degree (NOT ORDER):  ');   %polynomial degree

%%%%% UNCOMMENT THESE TWO LINES FOR UNIFORM INTERVALS
h = 1.0/N;
x = 0:h:1;
%%%%% UNCOMMENT THIS SINGLE LINE FOR NON_UNIFORM INTERVAL
% x = nonUniformIntervals( N, 1.0, 0 ); %element definitions

[z,w] = JacobiGZW(M+1,0.0,0.0);  %get quadrature points and weights
% six quadrature points at which to evaluate the funcitons -- non
% superconvergent
[zEval, wEval] = JacobiGZW(evalPoints, 0.0, 0.0);
NumP = size(zEval, 1);


for i = 0:M,
  Lvalue(i+1,:) = sqrt(i+0.5)*JacobiPoly(i,z,0.0,0.0)';
  LvalueMass(i+1) = dot(w,Lvalue(i+1,:).*Lvalue(i+1,:));
end
