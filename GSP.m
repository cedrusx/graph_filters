function [frq, Basis, Shift] = GSP(W, METHOD)
% W: weight matrix of a graph
% frq: array of frequencies
% Basis: basis matrix of graph Fourier transform
% Shift: shift matrix in vertex domain

if nargin < 2 || METHOD == 1 % 1: Laplacian (default)
    D = diag(sum(W,2));
    Shift = D - W;
elseif METHOD == 2           % 2: Normalized Laplacian
    D = diag(sum(W,2));
    Shift = eye(size(D)) - D^(-0.5)*W*D^(-0.5);
elseif METHOD == 3           % 3: Adjacent matrix
    Shift = W;
else
    error('Invalid method!');
end
[Basis,frq] = eigsort(Shift);
