function A = bestint(A,threshold)
% Make elements in A that are close to integers be integers
if nargin < 2
    threshold = 1e-10;
end
% A = A .* (abs(A) > threshold);
Aint = round(A);
mask = abs(A - Aint) < threshold;
A(mask) = Aint(mask);