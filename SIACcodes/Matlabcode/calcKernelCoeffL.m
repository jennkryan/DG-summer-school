function out = calcKernelCoeffL( deg, offset )
% calcCoeff calculates the kernel coefficients of degree with offset.
%   the degree is the degree of the BSPLINES not the polynomial against
%   which the kernel will be convolved.  Given degree, deg, there are
%   2*k-1 coeffcients which align from "left" to "right".  
%
%   The offset determines how the coefficients are centered.  0 centers
%   negative integers shifts the kernel one element to the left, and
%   positive values shift the kernel to the right

matSize = 2 * deg - 1;
k = deg - 1;
gamma_1 = -k + offset;
gamma_n = k + offset ;
A = zeros(matSize);

[z, w] = JacobiGZW( deg, 0.0, 0.0);
basis = zeros(matSize,1);
for i = 1:matSize
    basis(i) = JacobiPoly(i-1,0,0.0,0.0);
end

bs = getBSplinePP( deg );
kernelWidth = (3 * deg - 2);    % order + 2 * (order - 1)
halfWidth = kernelWidth * 0.5;
compactSupport = [-halfWidth, halfWidth] ;
kernelBreaks = (-halfWidth:halfWidth); 

% map legendre interval [-1, 1] to the full kernel width
legendreToKernel = kernelWidth / 2;  % 
kernelToLegendre = 1 / legendreToKernel;
[bs_breaks, c, d, l, k] = unmkpp(bs);

% precalculate legendre mapped to the Kernel domain
bp = 1;
%                   interval            quad points   L_i
yVals = zeros(size(kernelBreaks,2) - 1, size(z,1), matSize);
for r = 1:matSize
    power = r - 1;
    for i = 1:(size(kernelBreaks,2) - 1)
        a = kernelBreaks(i) * kernelToLegendre;
        b = kernelBreaks( i + 1) * kernelToLegendre;
        legendreZ = (z' + 1) * 0.5 * (b - a) + a;
        yVals(i,:, r) = JacobiPoly(power,legendreZ,0.0,0.0);
    end
end

    

for row = 1:matSize
    power = row - 1;
    gammaIndex = 1;
    for gamma = gamma_1:gamma_n
        % perform integral from -breakpoint to +breakpoint of
        intVal = 0;
        breaks = bs_breaks + gamma;
        for i = 1:(size(breaks,2) - 1)
            a = breaks( i );
            b = breaks( i + 1 );
            LIndex = gamma - gamma_1 + i;
            yVal = yVals(  LIndex, :, row);
            % mapped to kernel
            quadZ = (z' + 1) * 0.5 * (b - a) + a;
            bsVal = evalMappedPP( bs, 1.0, gamma, quadZ );
            %   perform quadrature
            quadVal = (b-a)*0.5*dot( w, bsVal .* yVal );
            intVal = intVal + quadVal;   
        end
        A(row, gammaIndex) = intVal;
        gammaIndex = gammaIndex + 1;
    end
end



out = (inv(A)*basis)';