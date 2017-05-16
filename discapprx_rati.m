function [nu,de] = discapprx_rati(x, y, order, METHOD, de_init)
% Polynomial approximant for a discrete function y = f(x)

nu = zeros(order(1) + 1, 1);
de = zeros(order(2) + 1, 1); de(1) = 1;
if order(1)+1 >= length(x)
    nu = discapprx_poly(x, y, order(1)); return;
elseif sum(order)+1 > length(x)
    order(2) = length(x) - order(1) - 1;
end

if ~exist('METHOD','var')
    METHOD = 2;
end

if METHOD == 1
% Method 1
% [nu,de] = argmin sum_square(y * polynomial(de,x) - polynomial(nu,x))
    A = ones(length(x), sum(order) + 1);
    for k = 1:order(1)
        A(:,k+1) = x(:).^k;
    end
    for k = 1:order(2)
        A(:,order(1)+1+k) = - x(:).^k .* y(:);
    end
    co = pinv(A,eps) * y(:);
    nu = co(1:(order(1)+1));
    de(2:(order(2)+1)) = co((order(1)+2):(order(1)+order(2)+1));
elseif METHOD == 2
% Method 2
% Carlos F. Borges, A Full-Newton Approach to Separable Nonlinear Least
% Squares Problems and its Application to Discrete Least Squares Rational
% Approximation, Electronic Transactions on Numerical Analysis, Volume 35,
% pp.57-68, 2009.
    if exist('de_init', 'var')
        [alpha, c] = dlsqrat(x(:),y(:),order(1),order(2),de_init(2:length(de_init)));
    else
        [alpha, c] = dlsqrat(x(:),y(:),order(1),order(2));
    end
    nu = c;
    de(2:(order(2)+1)) = alpha;
end