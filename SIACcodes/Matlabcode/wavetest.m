function [out] = wavetest(x,order,a,z,delta,alpha,beta)
% function [u] = wavetest(x,a,z,delta,alpha,beta)
% Purpose: Initial conditions for linear wave test problems.


N = length(x)-1; 
[zg,wg] = JacobiGZW(2*order,0.0,0.0);

LegPoly = zeros(order,2*order);
LegPoly(1,:) = 1;
if order >= 2
    LegPoly(2,:) = zg;
    for m=2:order-1
        LegPoly(m+1,:) = ((2*m-1).*zg'.*LegPoly(m,:)-(m-1.).*LegPoly(m-1,:))/(m);
    end
end
for m=1:order
    LegPoly(m,:) = sqrt(m-0.5).*LegPoly(m,:);
end
out = zeros(N,order);

for i=1:N 
    h = x( i + 1 ) - x( i );
    zmap = 0.5*h*(zg+1.0) + x(i);  %affine mapping of [-1,1] quadrature points on to element  
    [f] = fgauss(zmap,a,z,delta,alpha,beta);
    for mo = 0:order-1
        out(i,mo+1) =  dot(wg,f.*LegPoly(mo+1,:));
    end
end
return
% 