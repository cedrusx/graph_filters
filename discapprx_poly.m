function [co,preco] = discapprx_poly(x, y, order, METHOD)
% Polynomial approximant for a discrete function y = f(x)
% co = argmin sum_square(y - polynomial(co,x))

if ~exist('METHOD','var')
    METHOD = 1;
end

% Deduplication
[x,ix,~] = unique(x);
y = y(ix);

co = zeros(order + 1, 1);
if order >= length(x)
    order = length(x) - 1;
end

if METHOD == 0
% Method 0: MATLAB built-in function polyfit()
    co((order+1):-1:1) = polyfit(x, y, order);
elseif METHOD == 0.1
    [p,~,mu] = polyfit(x, y, order);
%     syms symy;
%     co((order+1):-1:1) = sym2poly(subs(poly2sym(p),(symy-mu(1))/mu(2)));
    co((order+1):-1:1) = p;
    preco = [-mu(1);1]/mu(2);    
elseif METHOD == 1
% Method 1: LS
% y = c0 + c1 * x + c2 * x^2 + ... + ck * x^k
    A = ones(length(x), order + 1);
    for k = 1:order
        A(:,k+1) = x(:).^k;
    end
%     co(1:(order+1)) = pinv(A,eps) * y(:);
% 上式仅在order等于11或12时略好于下式，大于12时远远差与下式，小于11时结果相同
    co(1:(order+1)) = A \ y(:);
elseif METHOD == 1.1
    A = fliplr(vander(x));
    A = A(:, 1:(order+1));
    co(1:(order+1)) = A \ y(:);
elseif METHOD == 1.2
    A = double(feval(symengine,'linalg::vandermonde',x));
    A = A(:, 1:(order+1));
    co(1:(order+1)) = pinv(A) * y(:);
elseif METHOD == 2
% Method 2
% y = c0 + c1(x-x1) + c2(x-x1)(x-x2) + ... + ck(x-x1)...(x-xk)
    A = ones(length(x), order + 1);
    for k = 1:order
        for kk = 1:k
            A(:,k+1) = A(:,k+1) .* (x(:) - x(kk));
        end
    end
    coa = pinv(A,eps) * y(:);
    syms s; fs = coa(1);
    for k = 1:order
        term = coa(k+1);
        for kk = 1:k
            term = term * (s - x(kk));
        end
        fs = fs + term;
    end
    co((order+1):(-1):1) = sym2poly(fs);
end