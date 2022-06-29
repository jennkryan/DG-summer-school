



% Calculate the exact solution for the linear transport equations
% put the characteristic in [-1,1]
xi = x-wavesp*FinalTime;
for i=1:length(xi(:))
    while ((xi(i) > xmax) || (xi(i)<xmin))
        if xi(i) > xmax
            xi(i) = xi(i)-L;
        elseif xi(i) < xmin
            xi(i) = xi(i)+L;
        end
    end
end

uetemp = uexact(xi(:));
ue = reshape(uetemp,m+1,N);

F = 14; % Font size

% Output for Linear Advection 
figure(1)
plot(x,ue,'k','LineWidth',2,'DisplayName','Exact')
hold on
plot(x,uh,':','LineWidth',3,'DisplayName','DG Approx')
xlabel('X','FontSize', F)
ylabel('u(x,T), u_h(x,T)','FontSize', F)
title(['Solution and Approximation for p=',num2str(m),...
    ' and N=',num2str(N),' elements'],'FontSize', 16)

figure(2)
plot(x,ue-uh,'LineWidth',3)
xlabel('X','FontSize', F)
ylabel('Error','FontSize', F)
title(['Error in approximation for p=',num2str(m),...
    ' and N=',num2str(N),' elements'],'FontSize', 16)

figure(3)
semilogy(x,abs(ue-uh),'LineWidth',3)
xlabel('X','FontSize', F)
ylabel('log(|Error|)','FontSize', F)
title(['Semi-log Errors for p=',num2str(m),...
    ' and N=',num2str(N),' elements'],'FontSize', 16)
