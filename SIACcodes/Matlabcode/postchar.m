function out = nonUniformPost( elements, boundaries, evalPoints, kernelScale )

% polynomials are all of degree k
%   map k+1 eval points (quadrature zeros) to each interval
%   that means there will be (k+1) * numElements = N points
%   output is a 2 X N array, first row contains x values in world coords
%   second row contains the corresponding y values

% set up criteria
[numElements, order] = size(elements);
approxorder = order;
kapprox = approxorder -1;
order = 0;
kpp =  order - 1;
[z, w] = JacobiGZW(evalPoints, 0.0, 0.0);

% set x values -- this could alternatively be supplied externally
sampleCount = numElements * evalPoints; %order;
out = zeros(2, sampleCount );
offset = 1;
last = offset + evalPoints-1;
for i = 1:numElements
    h = boundaries( i + 1) - boundaries(i);
    out(1, offset:last) = h * 0.5 * (z' + 1) + boundaries( i );
    offset = last + 1;
    last = offset + evalPoints-1;
end

% USE THIS IF WE WANT TO USE THE LARGEST ELEMENT IN THE MESH FOR THE KERNEL
% SCALING.
%kernelScale = 0;
if kernelScale ==0
    for n=1:numElements
        kernelScale = max(kernelScale,boundaries( n + 1 ) - boundaries( n ));
        elementWidth(n) = boundaries(n+1) - boundaries(n);
    end
end

% set up kernel statistics -- bspline for evaluation, etc.
bs = getBSplinePP( order );     % this is ORDER of the bspline -- i.e. order 1 = degree 0 (constant bspline)


%coeff = getKernelCoeffs( order );
%coeff = calcKernelCoeffL( order, 0 );
 % kernel width defined in KERNEL SPACE - where the width of bspline of degree 0 is 1
%kernelWidth = (3 * order - 2) ;   % order + 2 * (order - 1)
%halfWidth = kernelWidth * 0.5;
%compactSupport = [-halfWidth, halfWidth] ;
%kernelBreaks = (-halfWidth:halfWidth);

RS = max(ceil(0.5*(kapprox+order-1)),ceil(0.5*kapprox));
% convolution
kgl = ceil(0.5*(kapprox+order+1));
[z, w] = JacobiGZW(kgl, 0.0, 0.0);  % need enough quad points to support 2k + 1 degree
outIndex = 1;   % where to write the x value to in the output
kwidehalf = 0.5 * (2*RS+ell);
for n = 1:numElements
    elementWidth = boundaries( n + 1 ) - boundaries( n );
    % USE THIS IF WE WANT TO USE THE CURRENT ELEMENT SIZE FOR THE ELEMENT
     % SCALING.
%     kernelScale = elementWidth(n);
%     rtest = 4*k;
     %rtest = 2*k;
     %USE THIS PART FOR MAXIMUM SCALING WITHIN THE KERNEL SUPPORT
%     if (n <= rtest + order)  
%         kernelScale = 0;
%         for j=1:(rtest + order)
%             kernelScale = max(kernelScale,elementWidth(j));
%         end
%     elseif (numElements -n <= rtest + order)
%         kernelScale = 0;
%         for j=(numElements):-1:(numElements - (rtest+order))
%             kernelScale = max(kernelScale,elementWidth(j));
%         end
%     else
%         kernelScale = 0;
%         for j=(n-ceil(0.5*(3*k+1))):(n+ceil(0.5*(3*k+1)))
%             kernelScale = max(kernelScale,elementWidth(j));
%         end
%     end
    for xIndex = 1:evalPoints
        % find all the breakpoints for integration
        offset = out(1, (n - 1) * evalPoints + xIndex);  % get the breakpoint
        if (offset < boundaries(1) + 0.5 * (3 * order - 2)*kernelScale)
             % Near the left boundary -- use one-sided filtering
             rspline = 4*kpp;
             %rspline = 2*k;
             lambda = min(0,-0.5*(rspline+order)+(offset - boundaries(1))/kernelScale); 
         elseif (offset > boundaries(numElements+1) - 0.5 * (3 * order - 2)*kernelScale)
             % Near the right boundary -- use one-sided filtering
             rspline = 4*kpp;
             %rspline = 2*k;
             lambda = max(0,0.5*(rspline+order)+(offset - boundaries(numElements+1))/kernelScale);
         else
             % Domain interior -- symmetric post-processing
             rspline = 2*kpp;
             lambda = 0;        
         end
     %  Kernel stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
         coeff = KernelCoeffs( rspline, kpp , lambda);
         kernelWidth = rspline + order;
         halfWidth = kernelWidth * 0.5;
         compactSupport = [lambda - halfWidth, lambda + halfWidth] * kernelScale;
         kernelBreaks = (lambda - halfWidth:lambda + halfWidth) * kernelScale;

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        localCompactSupport = compactSupport * elementWidth + offset;   % scale compact support to element width
%        localKernelBreaks = kernelBreaks * elementWidth + offset;  % world coordinates
        if (offset >= boundaries(1) + 0.5 * (3 * order - 2)*kernelScale) && ...
                 (offset <= boundaries(numElements+1) - 0.5 * (3 * order - 2)*kernelScale)
             localCompactSupport = offset + compactSupport ;   % scale compact support to element width
             localKernelBreaks = offset + kernelBreaks ;  % world coordinates
        else
             localCompactSupport = offset - compactSupport ;   % scale compact support to element width
             localKernelBreaks = offset - kernelBreaks ;  % world coordinates
        end
        uhatBreaks = funcBreakPoints( boundaries, localCompactSupport );  % in world coordinates (outside of interval definition)
        breaks = [localKernelBreaks uhatBreaks];
        breaks = sort(breaks);

        % integrate product of convolved kernel and uhat
        intVal = 0;
        for i = 1:(size(breaks,2) - 1)
            a = breaks( i );
            b = breaks( i + 1 );
            %   map quadrature points
            quadZ = (z' + 1) * 0.5 * (b - a) + a;
            % evaluate uhat and kernel at quadrature points
            uhatVal = periodicEval( elements, boundaries, quadZ );
            kernelVal = evalKernelPP( bs, rspline, kpp, lambda, elementWidth, offset, coeff, quadZ);
            %   perform quadrature
            intVal = intVal + (b-a)*0.5*dot( w, kernelVal .* uhatVal );   
        end
        out(2,outIndex) = intVal / elementWidth;  % rescale the integral by the element width
        outIndex = outIndex + 1;
    end
end
    