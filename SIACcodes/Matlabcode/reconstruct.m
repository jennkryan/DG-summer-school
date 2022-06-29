function urecon = reconstruct(b,xq)
%Reconstruct u based on b

[Nmesh modes] = size(b);
%Lq = Legendre_Polynomials(k,xq);
for i = 1:modes
  Lq(i,:) = sqrt(i-0.5)*JacobiPoly(i-1,xq,0.0,0.0)';
end

utemp=[];
qpts = length(xq);
for j=1:Nmesh
    for j2 = 1:qpts
        uapprox = 0;
        for ell=1:modes
            uapprox = uapprox+b(j,ell)*Lq(ell,j2); % DG approximation in points xq (local coordinates)
        end
        utemp = [utemp uapprox];
    end
end
urecon = utemp;
