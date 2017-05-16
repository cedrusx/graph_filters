function y = rational(nu,de,x)

if min(size(x)) == 1
    y = polynomial(nu,x) ./ polynomial(de,x);
else
    y = polynomial(de,x) \ polynomial(nu,x);
end