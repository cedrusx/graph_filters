function c = filterco(f,lambset)
% f: filter coefficients, s.t. F(y) = F(x) .* f
% lambset: N-by-1 array or N-by-N diagonal matrix containing frequencies
% PURPOSE: y = polynomial(c,L) * x
N = length(f);
if min(size(lambset)) > 1
    lambset = diag(lambset);
end
lambset = reshape(lambset,[],1);
A = zeros(N);
for k = 1:length(f)
    A(:,k) = lambset.^(k-1);
end
c = A\f;
