function A=matLattice(n,loop)
% ����n��n������������Ӿ���
% loop=1��ʾ�����ޱ߽磨���ұ����������±�������

nn=n*n;
A=zeros(nn,nn);
for k=1:nn
    if mod(k,n)~=0
        A(k,k+1)=1;A(k+1,k)=1;
    end
end
for k=1:(nn-n)
    A(k,k+n)=1;A(k+n,k)=1;
end

if loop
    for k=n:n:nn
        A(k,k-n+1)=1;A(k-n+1,k)=1;
    end
    for k=1:n
        A(k,k+nn-n)=1;A(k+nn-n,k)=1;
    end
end