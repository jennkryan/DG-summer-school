function [out] = initprojNonUni(order,Elmts,x,func)

[z,w] = JacobiGZW(2*order,0.0,0.0);

ONLegPoly = zeros(order,2*order);
for i=1:order
    ONLegPoly(i,:) = sqrt(i-0.5)*legendreP(i-1,z);
end


out = zeros(Elmts,order);
for nel = 1:Elmts
        h = x( nel + 1 ) - x( nel );
        zmap = 0.5*h*(z+1.0) + x(nel);  %affine mapping of [-1,1] quadrature points on to element
        f = feval( func, zmap);           %function to be mapped evaluated at mapped quadrature points
        for i = 0:order-1
            out(nel,i+1) =  dot(w,f'.*ONLegPoly(i+1,:));
        end
    end

    save('initialmodes.txt','out','-ascii')
end
