function y = polynomial(c,x)
% For scalar or array:
% y = c(1) + c(2)*x + c(3)*x.^2 + ...
% For square matrix:
% y = c(1)*I + c(2)*x + c(3)*x^2 + ...
% Another choice: polyval()

y = zeros(size(x));
if min(size(x)) == 1 % scalar or array
    for k = 1:length(c)
        y = y + x.^(k-1) * c(k);
    end
elseif min(size(x)) == max(size(x)) % square matrix
    for k = 1:length(c)
        y = y + x^(k-1) * c(k);
    end
else
    fprintf('Error using polynomial()\n');
end