function [snu,sde,tail] = partialfraction(nu,de,domain)
% domain == 'R': Partial faction composition over the real set (default)
% domain == 'C': Partial faction composition over the complex set

if nargin > 2 && domain == 'C'
    [res,poles,tail] = residue(nu(end:-1:1),de(end:-1:1));
    snu = cell(size(res));
    sde = cell(size(res));
    for k = 1:length(res)
        sde{k} = [1, -1/poles(k)];
        snu{k} = [-res(k)/poles(k)];
    end
    tail = tail(end:-1:1);
    return
end

syms x;
nu = nu(:)'; de = de(:)';
y = sum(nu.*x.^((1:length(nu))-1))/sum(de.*x.^((1:length(de))-1));
y = feval(symengine,'partfrac', y, x, 'Domain = R_');
ys = children(y);
if ~isnan(0/(sum(ys)-y)) 
    % 判断ys是否是根据加法分解的
    % 若不是，说明y只有一项，另ys=y
    ys = y;
end
snu = cell(0);
sde = cell(0);
tail = [];
for k = 1:length(ys)
    [n,d] = numden(ys(k));
    [sde1,deg] = mysym2poly(d,x);
    snu1 = mysym2poly(n,x);
    if deg > 0 % 分母阶数大于0，则该项作为分解项
        sde = [sde;sde1];
        snu = [snu;snu1];
    else % 分母阶数等于0，则该项并入tail
        if length(snu1) > length(tail)
            tail = [tail,zeros(1,length(snu1)-length(tail))];
        end
        tail(1:length(snu1)) = tail(1:length(snu1)) + snu1/sde1;
    end
end
return

function [co,deg] = mysym2poly(f,x)
deg = feval(symengine,'degree',f,x);
if deg > 0  % 判断阶数是因为当f只有常数项时不能用sym2poly函数
    co = sym2poly(f); 
    co = co(end:-1:1);
else
    co = [eval(f)];
end