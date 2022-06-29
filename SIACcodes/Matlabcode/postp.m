%ell=M+1; % order of B-spline for filtering  For superconvergence should be M+1
ellp2=input('Input smoothness required (>=0).  0 = continuous:  ');
ell = ellp2+2;
%!!!!!!!!!!!!!!! IMPORTANT!!!!!!!!!!!!!!!!!!!!!!


tic;        % time it
% set kernelScale = 0 in order to use the largest element.
kernelscale = input('Input Kernel Scaling: 0 = element width; 0.5 for half-width; 2 for double width:  ');
kernelScale = kernelscale*h;
if ellp2 >= 0
    filtered = nonUniformPost(uhat, x, evalPoints, ell, kernelScale );
elseif ellp2 == -1
    filtered = postchar(uhat, x, evalPoints, kernelScale);
else
    fprintf('ERROR:  Negative B-spline order.')
end
toc;