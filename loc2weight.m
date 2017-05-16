function w = loc2weight(Loc, A, theta)

if ~exist('theta', 'var')
    theta = 1;
end
dis = loc2dis(Loc, A);
w = exp(-(dis/theta).^2/2);

w = w - diag(diag(w));