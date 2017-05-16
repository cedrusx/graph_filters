function [U,evalue] = eigsort(A)
[U,evalue] = eig(A);
evalue = diag(evalue);
[evalue,ord] = sort(evalue,'ascend');
U = U(:,ord);