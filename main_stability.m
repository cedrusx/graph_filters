% A = matLocalRConnected(10,0.3);
% [frq, U, L] = GSP(A,1);
frq = [0.0000    0.2496    1.1572    1.9941    3.0000    ...
    3.4254    4.0000    4.4336    5.4560    5.5841]';
cutoff = 3;

h1 = @(f)1./(1+(f/cutoff).^4);
% hi = @(x)(x<cutoff)+0;
% [nu, de] = discapprx_rati(0:maxfrq, hi(0:maxfrq), [0,3], 2);
% h2 = @(f)rational(nu,de,f);
h2 = @(f)1./(1+2*f-1*f.^2);

xlim = [0,2*cutoff];
ylim = [-1,1.5];
x = 0:0.01:xlim(2);
y1 = h1(x); y2 = h2(x);
[~,pole] = max(abs(y2));
pole = x(pole);
y2(abs(y2)>10) = NaN;

figure; hold on;
plot(xlim,[0,0],'-k');
hfilter = plot(x,y1,'-.r',x,y2,'-b','LineWidth',2);
plot([pole,pole],ylim,':b');
hpole = plot(pole,0,'xb','MarkerSize',10);
% stem(frq,h1(frq),':','Marker','None');
% stem(frq,h2(frq),':','Marker','None');
hfrq = plot(frq,zeros(size(frq)),'ok');
% hfig = plot(frq(1),h1(frq(1)),'-o',frq(1),h2(frq(1)),'-s',frq(1),h1(frq(1)),':','LineWidth',2);
legend([hfilter;hpole;hfrq], 'Filter 1', 'Filter 2', 'Poles', 'Spectrum components');
xlabel('Spectrum');
ylabel('Spectral response');
axis([0,max(frq),ylim]);