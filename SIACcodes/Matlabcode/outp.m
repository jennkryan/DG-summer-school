%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOT TRUE DATA AND APPROXIMATION DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



count = 1;
for nel = 1:N,
  h = x( nel + 1 ) - x( nel );
  zmap = 0.5*h*(zEval'+1.0) + x(nel);
  f = func(zmap);            %function to be mapped evaluated at mapped plotting points
  
  xvec(count:count+NumP-1) = zmap;
  ftrue(count:count+NumP-1) = f;
  
  %reconstruct polynomial as linear combination of basis functions
  utemp = zeros(size(zEval));
  for i = 0:M
    utemp = utemp + uhat(nel,i+1)*sqrt(i+0.5)*JacobiPoly(i, zEval, 0.0 ,0.0); %Lvalue(i+1,:)';
  end
  fapprox(count:count+NumP-1) = utemp;
  count = count + NumP;
end

save('fapprox.txt','fapprox','-ascii');
% plot 
error = max(abs(ftrue-fapprox));
close all;
figure(1)
plot(xvec,ftrue,'-b','LineWidth',2); hold on; grid on
plot(xvec,fapprox,'r--','LineWidth',2);hold on;
intPoints = feval( func, x );
plot(x, intPoints, 'bo');
