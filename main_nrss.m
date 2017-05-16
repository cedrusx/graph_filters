clear;
N = 1000; radius = 0.06; maxfrq = 25;
ord = 4;
cutoff = 10;
tmax = 50;
nu = [1];
de = [1 zeros(1,ord-1) 1/cutoff.^ord];
step = 2/sum(polynomial(de,[0 maxfrq]));
repmax = 2;
nrss_fast0 = zeros(repmax,tmax);
nrss_fast1 = zeros(repmax,tmax);
nrss_fast2 = zeros(repmax,tmax);
bound_fast = zeros(repmax,tmax);
tic;
for repid = 1:repmax
    while 1
        A = matLocalRConnected(N,radius);
        [frq, U, L] = GSP(A,1);
        if max(frq) < maxfrq
            break
        end
    end
    x = randn(N,1);
    y = IIR(nu,de,L,x,0);
%     y = U * ((frq<cutoff) .* (U' * x));
    nrss = @(x)(y-x)'*(y-x)/(y'*y);
    gfrq = polynomial(de,frq);
    for t = 1:tmax
        y0 = IIR(nu,de,L,x,2,t,0,step);
        nrss_fast0(repid,t) = nrss(y0);
        y1 = IIR(nu,de,L,x,2,t,1,0,[0 maxfrq]);
        nrss_fast1(repid,t) = nrss(y1);
        y2 = IIR(nu,de,L,x,2,t,2,0,[0 maxfrq]);
        nrss_fast2(repid,t) = nrss(y2);
        bound_fast(repid,t) = max((gfrq-1).^2.*(1-step*gfrq).^(2*t));
    end
end
toc;
figure;
comm = (1:tmax);
semilogy(comm,mean(nrss_fast0),'-.r',...
    comm,mean(nrss_fast1),'--b',comm,mean(nrss_fast2),'-k',...
    comm,mean(bound_fast),':r','LineWidth',1.5);
legend('Direct form','Cascade form','Parallel form',...
    'Upper bound of direct form');
xlabel('Number of iterations');
ylabel('NRSS');
grid on;