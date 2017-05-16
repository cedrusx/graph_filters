function y = ChebFilter(h,L,x,M,range)
% Distributed graph filtering via Chebyshev polynomial approximation
% h: filter's spectral response (function handle)
% L: Laplacian matrix
% x: input signal 
% M: maximal order of Chebyshev polynomials
% range = [minfrq maxfrq]: bounds on the spectrum of the graph
% y: filtered signal

if ~exist('range','var')
    range = [min(eig(L)),max(eig(L))];
    % In the script "errcompfull", L has been shifted before filtering
end

% Step 1. Calculation of coefficients
c = chebcoeffs(chebfun(h,range,M));

% Step 2. Distributed filtering
y = chebpolyvalm(c,L,range) * x;
