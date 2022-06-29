%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOT TRUE DATA AND APPROXIMATION DATA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fapproxtemp = reconstruct(uhat,zEval);
fapprox = reshape(fapproxtemp,evalPoints,N);
%count = 1;
%for nel = 1:N
%  h = x( nel + 1 ) - x( nel );
%  zmap = 0.5*h*(zEval'+1.0) + x(nel);
%  f = func(zmap);            %function evaluated at mapped plotting points
  
%  xvec(count:count+NumP-1) = zmap;
%  ftrue(count:count+NumP-1) = f;
  
%  %reconstruct polynomial as linear combination of basis functions
%  utemp = zeros(size(zEval));
%  for i = 0:M
%    utemp = utemp + uhat(nel,i+1)*JacobiPoly(i, zEval, 0.0 ,0.0); %Lvalue(i+1,:)';
%  end
%  fapprox(count:count+NumP-1) = utemp;
%  count = count + NumP;
%end

%save('fapprox.txt','fapprox','-ascii');
% plot 
error = max(abs(fmapEval-fapprox));
close all;
figure(1)
plot(zmapEval,fmapEval,'-b','LineWidth',2); hold on; grid on
plot(zmapEval,fapprox,'r--','LineWidth',2);hold on;
intPoints = func(x);
plot(x, intPoints, 'bo');
