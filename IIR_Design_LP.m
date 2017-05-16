function [nu,de,order,fc] = IIR_Design_LP(dp,ds,fp,fs)
% Analytical IIR low-pass filter design 
% h(f) = 1./(1+(f/f0).^order) = polynomial(f,nu)./polynomial(f,de)
% Specification:
%  1-dp <= h(f) <= 1+dp, 0 <= f <= fp
%  -ds <= h(f) <= ds, f >= fs

order = ceil(log((1/(1-dp)-1)/(1/ds-1)) / log(fp/fs));
fc = fp/(1/(1-dp)-1)^(1/order);
nu = 1;
de = [1 zeros(1,order-1) 1/fc.^order];