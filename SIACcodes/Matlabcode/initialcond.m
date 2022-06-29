%func = @(x) ones(size(x));
%func = @(x) 1-abs(x));
%func = @(x) 1-x.^2;
func = @(x) sin(4*pi*x);


fmapquad = func(zmapquad);
fmapEval = func(zmapEval);

%calculate modes
for nel = 1:N
  for i = 0:M
      quadsum = 0;
      for ni = 1:M+1
          quadsum = quadsum + w(ni)*fmapquad(ni,nel).*Lvalue(i+1,ni);
      end
    uhat(nel,i+1) =  quadsum;
  end
end

