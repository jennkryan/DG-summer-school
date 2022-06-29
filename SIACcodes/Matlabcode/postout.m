figure(1)
plot(filtered(1,:), filtered(2,:), '-go','LineWidth',2);
legend('true','approx', 'interval markers', 'post');
s = sprintf('Discrete Linf Error = %0.5g',error);
xlabel(s);
title('Function, approximation and post processed');
% 
fapproxError = abs(ftrue-fapprox)+ 1.0e-14 ;
ftruepp = func(filtered(1,:));
filterError = abs(ftruepp - filtered(2,:))+ 1.0e-14 ;
figure(2);
title('Error- log 10');
semilogy(xvec, fapproxError,'r','LineWidth',2);
hold on;
semilogy(filtered(1,:), filterError, 'b','LineWidth',2);
legend('Approx error', 'Filter error');
ylabel('log_1_0( error )');
xlabel('x, world space');
title('Error of approximation and post process');

ppapprox = filtered(2,:);

s1 = sprintf('p=%2i,   ell = %2i\n',M,ell);
s2 = sprintf('      N    DG L-Inf Error    PostP DG L-inf Error\n');
s3 = sprintf('    %3i   %13.3d   %13.3d\n',N,max(abs(fapproxError)),max(abs(filterError)));
disp([s1,s2,s3])