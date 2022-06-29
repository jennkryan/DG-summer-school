%func = @one;
%func = @linear;
%func = @quadratic;
%func = @cubic;
func = @(x) sin(4*pi*x);

for nel = 1:N
  h = x( nel + 1 ) - x( nel );
  zmap = 0.5*h*(z+1.0) + x(nel);  %affine mapping of [-1,1] quadrature points on to element
  f = func(zmap);           %function to be mapped evaluated at mapped quadrature points
  for i = 0:M
    uhat(nel,i+1) =  dot(w,f'.*Lvalue(i+1,:))/LvalueMass(i+1);
  end
end

