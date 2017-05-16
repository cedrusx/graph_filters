% 2014.12.6 for SPL revision
% With disturbance of network topology
% With distributed realization of IIR
% Spectrum is normalized in approximation
clear;
DoCheb = 1; % Chebyshev
DoIDIIR = 0; tm = 100; % IDIIR with tm iterations
DoCDIIR = 0; M = 10; % IIR by M-degree Chebyshev approximation
DoTmax = 0; kappa = 0.01; % tmax (direct form)

N = 100; radius = 0.2;
ordset = 1:17;
cutoff = 10;
h = @(x)((x<cutoff)+0); % Low pass
repmax = 100;

rmse_poly1 = zeros(repmax,length(ordset));
rmse_poly2 = zeros(repmax,length(ordset));
rmse_rati21 = zeros(repmax,length(ordset));
rmse_rati22 = zeros(repmax,length(ordset));
rmse_rati41 = zeros(repmax,length(ordset));
rmse_rati42 = zeros(repmax,length(ordset));
if DoCheb
rmse_cheb1 = zeros(repmax,length(ordset));
rmse_cheb2 = zeros(repmax,length(ordset));
end
if DoIDIIR
rmse_rati21d = zeros(repmax,length(ordset));
rmse_rati22d = zeros(repmax,length(ordset));
rmse_rati41d = zeros(repmax,length(ordset));
rmse_rati42d = zeros(repmax,length(ordset));
end
if DoTmax
tmax2 = zeros(repmax,length(ordset));
tmax4 = zeros(repmax,length(ordset));
end
if DoCDIIR
rmse_rati21c = zeros(repmax,length(ordset));
rmse_rati22c = zeros(repmax,length(ordset));
rmse_rati41c = zeros(repmax,length(ordset));
rmse_rati42c = zeros(repmax,length(ordset));
end
tic;
parfor repid = 1:repmax
    [A1, Loc] = matLocalRConnected(N,radius);
    [frq1, U1, L1] = GSP(A1,1);
    hf1 = h(frq1);
    A2 = modifynetwork(A1,0.1,3,Loc,radius);
    [frq2, U2, L2] = GSP(A2,1);
    x = randn(N,1);
    y1 = U1 * (h(frq1) .* (U1' * x));
    y2 = U2 * (h(frq2) .* (U2' * x));
    rmse_poly1s = zeros(1,length(ordset));
    rmse_poly2s = zeros(1,length(ordset));
    rmse_rati21s = zeros(1,length(ordset));
    rmse_rati22s = zeros(1,length(ordset));
    rmse_rati41s = zeros(1,length(ordset));
    rmse_rati42s = zeros(1,length(ordset));
    if DoIDIIR
    rmse_rati21ds = zeros(1,length(ordset));
    rmse_rati22ds = zeros(1,length(ordset));
    rmse_rati41ds = zeros(1,length(ordset));
    rmse_rati42ds = zeros(1,length(ordset));
    end
    if DoCDIIR
    rmse_rati21cs = zeros(1,length(ordset));
    rmse_rati22cs = zeros(1,length(ordset));
    rmse_rati41cs = zeros(1,length(ordset));
    rmse_rati42cs = zeros(1,length(ordset));
    end
    if DoTmax
    tmax2s = zeros(1,length(ordset));
    tmax4s = zeros(1,length(ordset));
    end
    if DoCheb
    rmse_cheb1s = zeros(1,length(ordset));
    rmse_cheb2s = zeros(1,length(ordset));
    end
    for ordid = 1:length(ordset)
        order = ordset(ordid);
        % FIR
        [co_poly,co_pre] = discapprx_poly(frq1, hf1, order, 0.1);
        fp = @(x)polynomial(co_pre,x);
        rmse_poly1s(ordid) = rmse(y1,polynomial(co_poly,fp(L1))*x);
        rmse_poly2s(ordid) = rmse(y2,polynomial(co_poly,fp(L2))*x);
        % Chebyshev
        if DoCheb
        rmse_cheb1s(ordid) = rmse(y1,ChebFilter(h,L1,x,order,[0 max(frq1)]));
        rmse_cheb2s(ordid) = rmse(y2,ChebFilter(h,L2,x,order,[0 max(frq1)]));
        end
        % IIR (denominator degree 2)
%         fp = @(x)x;
        [nu_rati, de_rati] = discapprx_rati(fp(frq1), hf1, [order,1], 2, [1; 0]);
        [nu_rati, de_rati] = discapprx_rati(fp(frq1), hf1, [order,2], 2, [de_rati; 0]);
        rmse_rati21s(ordid) = rmse(y1,IIR(nu_rati,de_rati,fp(L1),x));
        rmse_rati22s(ordid) = rmse(y2,IIR(nu_rati,de_rati,fp(L2),x));
        % IDIIR (denominator degree 2)
        if DoTmax
        tmax2s(ordid) = tmax(kappa,polynomial(de_rati,fp(frq1)),2);
        end
        if DoIDIIR
        rmse_rati21ds(ordid) = rmse(y1,IIR(nu_rati,de_rati,fp(L1),x,1,tm,2,0,fp(frq1)));
        rmse_rati22ds(ordid) = rmse(y2,IIR(nu_rati,de_rati,fp(L2),x,1,tm,2,0,fp(frq1)));
        end
        % Discrete IIR via Chebyshev (denominator degree 2)
        if DoCDIIR
        rmse_rati21cs(ordid) = rmse(y1,IIR(nu_rati,de_rati,fp(L1),x,3,M,0,0,fp(frq1)));
        rmse_rati22cs(ordid) = rmse(y2,IIR(nu_rati,de_rati,fp(L2),x,3,M,0,0,fp(frq1)));
        end
        % IIR (denominator degree 4)
        [nu_rati, de_rati] = discapprx_rati(fp(frq1), hf1, [order,3], 2, [de_rati; 0]);
        [nu_rati, de_rati] = discapprx_rati(fp(frq1), hf1, [order,4], 2, [de_rati; 0]);
        rmse_rati41s(ordid) = rmse(y1,IIR(nu_rati,de_rati,fp(L1),x));
        rmse_rati42s(ordid) = rmse(y2,IIR(nu_rati,de_rati,fp(L2),x));
        % IDIIR (denominator degree 4)
        if DoTmax
        tmax4s(ordid) = tmax(kappa,polynomial(de_rati,fp(frq1)),2);
        end
        if DoIDIIR
        rmse_rati41ds(ordid) = rmse(y1,IIR(nu_rati,de_rati,fp(L1),x,1,tm,2,0,fp(frq1)));
        rmse_rati42ds(ordid) = rmse(y2,IIR(nu_rati,de_rati,fp(L2),x,1,tm,2,0,fp(frq1)));
        end
        % Discrete IIR via Chebyshev (denominator degree 4)
        if DoCDIIR
        rmse_rati41cs(ordid) = rmse(y1,IIR(nu_rati,de_rati,fp(L1),x,3,M,0,0,fp(frq1)));
        rmse_rati42cs(ordid) = rmse(y2,IIR(nu_rati,de_rati,fp(L2),x,3,M,0,0,fp(frq1)));
        end
    end
    rmse_poly1(repid,:) = rmse_poly1s;
    rmse_poly2(repid,:) = rmse_poly2s;
    rmse_rati21(repid,:) = rmse_rati21s;
    rmse_rati22(repid,:) = rmse_rati22s;
    rmse_rati41(repid,:) = rmse_rati41s;
    rmse_rati42(repid,:) = rmse_rati42s;
    if DoIDIIR
    rmse_rati21d(repid,:) = rmse_rati21ds;
    rmse_rati22d(repid,:) = rmse_rati22ds;
    rmse_rati41d(repid,:) = rmse_rati41ds;
    rmse_rati42d(repid,:) = rmse_rati42ds;
    end
    if DoCDIIR
    rmse_rati21c(repid,:) = rmse_rati21cs;
    rmse_rati22c(repid,:) = rmse_rati22cs;
    rmse_rati41c(repid,:) = rmse_rati41cs;
    rmse_rati42c(repid,:) = rmse_rati42cs;
    end
    if DoTmax
    tmax2(repid,:) = tmax2s;
    tmax4(repid,:) = tmax4s;
    end
    if DoCheb
    rmse_cheb1(repid,:) = rmse_cheb1s;
    rmse_cheb2(repid,:) = rmse_cheb2s;
    end
    fprintf('Instance #%d finished.\n',repid);
end
toc;

figure;
ordp = 1:17;
pf = @(y,linespec)semilogy(ordset(ordp),mymean(y(:,ordp),1000),...
    linespec,'LineWidth',1.5);
ls = cell(0);
% Original network
if DoCheb
pf(rmse_cheb1,'-*m'); hold on; ls = [ls;'Chebyshev FIR'];
end
pf(rmse_poly1,'-^k'); ls = [ls;'Least-squares FIR'];
pf(rmse_rati21,'-ob'); ls = [ls;'IIR with denominator degree 2'];
pf(rmse_rati41,'-sr'); ls = [ls;'IIR with denominator degree 4'];
if DoIDIIR
pf(rmse_rati21d,'-og'); ls = [ls;['IDIIR with denominator degree 2 by ',num2str(tm),' iterations']];
pf(rmse_rati41d,'-sc'); ls = [ls;['IDIIR with denominator degree 4 by ',num2str(tm),' iterations']];
end
if DoCDIIR
pf(rmse_rati21c,'-og'); ls = [ls;['IIR with denominator degree 2 by ',num2str(M),'-degree Chebyshev']];
pf(rmse_rati41c,'-sc'); ls = [ls;['IIR with denominator degree 2 by ',num2str(M),'-degree Chebyshev']];
end
legend(ls);
% Perturbed network
if DoCheb
pf(rmse_cheb2,':*m');
end
pf(rmse_poly2,':^k');
pf(rmse_rati22,':ob'); 
pf(rmse_rati42,':sr'); 
if DoIDIIR
pf(rmse_rati22d,':og');
pf(rmse_rati42d,':sc');
end
if DoCDIIR
pf(rmse_rati22c,':og');
pf(rmse_rati42c,':sc');
end
xlabel('Degree of numerator'); ylabel('RMSE'); grid on;
% legend('FIR in original network',...
%     'IIR with denominator degree 2 in original network',...
%     'IIR with denominator degree 4 in original network',...
%     'FIR in perturbed network',...
%     'IIR with denominator degree 2 in perturbed network',...
%     'IIR with denominator degree 4 in perturbed network');

% figure;
% pf = @(x,y,linespec)semilogy(x,y,linespec,'LineWidth',1.5);
% pf(mymean(rmse_poly1),mymean(rmse_poly2),'^k'); hold on;
% pf(mymean(rmse_rati21),mymean(rmse_rati22),'ob'); hold on;
% pf(mymean(rmse_rati41),mymean(rmse_rati42),'sr'); hold on;
% xlabel('RMSE in original network');
% ylabel('RMSE in perturbed network');

if DoTmax
figure; 
ordp = 1:17;
pf = @(y,linespec)semilogy(ordset(ordp),mymean(y(:,ordp)),...
    linespec,'LineWidth',1.5);
pf(tmax2,'-og'); hold on;
pf(tmax4,'-sc');
legend('IIR with denominator degree 2','IIR with denominator degree 4');
xlabel('Degree of numerator'); ylabel('Required number of iterations');

figure; 
ordp = 1:17;
pf = @(y,linespec)semilogy(ordset(ordp),mymean(y(:,ordp)),...
    linespec,'LineWidth',1.5);
pf(tmax2*2,'-og'); hold on;
pf(tmax4*4,'-sc');
legend('IIR with denominator degree 2','IIR with denominator degree 4');
xlabel('Degree of numerator'); ylabel('Communication cost');
end